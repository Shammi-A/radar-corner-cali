function cal = corner_calibrate(cube, p, varargin)
% CORNER_CALIBRATE  Corner-reflector based channel calibration (AWR2243 + DCA1000)
% Inputs
%   cube : [numRx x Ns x Nc_pf x Nf] complex/real ADC (from dca1000_read_bin)
%   p    : RadarParams from mmws_parse_log (must include fields used below)
% Name-Value (common)
%   'cornerRange_m'   : known corner distance (meters) (recommended)
%   'rangeSearch_m'   : [rMin rMax] search window (m). If not given, use ±20 bins around estimate
%   'useFrames'       : number of frames to use (default: min(8, size(cube,4)))
%   'plot'            : true/false quick sanity plots (default: false)
%   'NfftR'           : range FFT size (default: nextpow2(2*Ns))
%   'NfftD'           : Doppler FFT size per-TX (default: nextpow2(chirpsPerLoop))
%
% Output struct 'cal' contains:
%   cal.rxWeights    : [numRx x 1] complex weights (apply per-RX)
%   cal.txWeights    : [numTxActive x 1] complex weights (if TDM; else empty)
%   cal.rangeBias_m  : scalar (meas - known)
%   cal.rangeBin     : detected range bin
%   cal.dR_m         : meters per range bin (with chosen NfftR)
%   cal.kdopp        : detected Doppler bin(s) (per-TX if TDM; else scalar)
%   cal.meta         : details incl. mapping and parameters

    ip = inputParser;
    ip.addParameter('cornerRange_m', []);
    ip.addParameter('rangeSearch_m', []);
    ip.addParameter('useFrames', []);
    ip.addParameter('plot', false);
    ip.addParameter('NfftR', []);
    ip.addParameter('NfftD', []);
    ip.parse(varargin{:});
    o = ip.Results;

    % ---- sizes & defaults ----
    [numRx, Ns, Nc_pf, Nf] = size(cube);
    if isempty(o.useFrames), o.useFrames = min(8, Nf); end
    useF = min(o.useFrames, Nf);
    cube = cube(:,:,:,1:useF);

    if isempty(o.NfftR)
        o.NfftR = 2^nextpow2( max(512, 2*Ns) );
    end
    if isempty(o.NfftD)
        o.NfftD = 2^nextpow2( max(32, max(1, p.chirpsPerLoop)) );
    end

    % ---- constants derived from config ----
    Fs     = p.Fs_Hz;
    S      = p.Slope_Hz_per_s;
    dR     = 3e8 * (Fs/o.NfftR) / (2*S);     % meters per range bin with NfftR

    % estimated bin from known corner range (if provided)
    if ~isempty(o.cornerRange_m)
        kEst = round(o.cornerRange_m / dR);
    else
        kEst = []; % we will scan for max
    end

    % ---- build chirp-to-TX mapping (TDM or not) ----
    [txOfChirp, patternLen, txActive, isTDM] = tx_mapping_from_params(p);
    numTxAct = numel(txActive);

    % ---- Range FFT for all chirps/frames (vectorized; preallocate to NfftR) ----
    % reshape to [numRx x Ns x (Nc_pf*useF)]
    X = reshape(cube, numRx, Ns, Nc_pf*useF);
    % window along fast time
    wR = hann(Ns,'periodic').'; wR = wR / max(eps, norm(wR));
    Xw = bsxfun(@times, X, reshape(wR, 1, Ns, 1));
    % range FFT along dim-2 with size NfftR
    Xr = fft(Xw, o.NfftR, 2);                         % [numRx x NfftR x Nc_pf*useF]

    % ---- choose range bin k0 ----
    if isempty(o.rangeSearch_m)
        if isempty(kEst)
            rSel = 1:floor(o.NfftR/2);               % search whole positive half
        else
            win = max(5, 20);                         % ±20 bins default
            rSel = max(1,kEst-win) : min(floor(o.NfftR/2), kEst+win);
        end
    else
        kMin = max(1, floor(o.rangeSearch_m(1)/dR));
        kMax = min(floor(o.NfftR/2), ceil(o.rangeSearch_m(2)/dR));
        if kMax < kMin, tmp=kMin; kMin=kMax; kMax=tmp; end
        rSel = kMin:kMax;
    end
    % energy across RX & chirps/frames on selected bins
    P_r = squeeze(sum(sum(abs(Xr(:,rSel,:)).^2,1),3)).';          % [numel(rSel) x 1]
    [~, idxMax] = max(P_r);
    k0 = rSel(idxMax);

    % ---- per-TX slow-time sequences at k0 (for stationary target) ----
    % slice range bin k0:  [numRx x 1 x Nc_pf*useF] -> [numRx x Nc_pf x useF]
    Xk = reshape( Xr(:,k0,:), numRx, Nc_pf, useF );

    % group chirps per-TX (TDM) OR use all chirps together (non-TDM)
    Htxrx = complex(zeros(max(1,numTxAct), numRx));   % accumulator per tx,rx
    kdopp = zeros(max(1,numTxAct),1);

    if isTDM
        % For each TX, pick its chirps across the frame (pattern repeats 'loops' times)
        loops = max(1, round(p.chirpsPerLoop));
        for tIdx = 1:numTxAct
            tx = txActive(tIdx);                                  % 1..3
            idxInPat = find(txOfChirp(1:patternLen) == tx);       % chirps in pattern for this TX
            if isempty(idxInPat), continue; end
            % collect chirps for this TX across the frame (pattern repeats)
            chirpSel = [];
            for ii = 1:numel(idxInPat)
                c0 = idxInPat(ii);
                chirpSel = [chirpSel, c0:patternLen:(c0+(loops-1)*patternLen)]; %#ok<AGROW>
            end
            chirpSel = chirpSel(chirpSel>=1 & chirpSel<=Nc_pf);

            % average across frames (stationary → DC), then coherent sum over chirps
            Xsel = mean( Xk(:,chirpSel,:), 3 );                   % [numRx x Nchirps(≈loops)]
            h = sum(Xsel, 2);                                     % [numRx x 1]
            Htxrx(tIdx,:) = h(:).';
            kdopp(tIdx)   = 0;                                    % DC used
        end
    else
        % Not TDM: cannot separate TX. Calibrate RX only using all chirps
        Xsel = mean( Xk, 3 );                                     % [numRx x Nc_pf]
        h = sum(Xsel, 2);                                         % DC across slow-time
        Htxrx(1,:) = h(:).';
        kdopp(1)   = 0;
    end

    % ---- RX weights (phase+mag) from boresight point target ----
    % reference = first enabled RX (LSB in mask)
    rxList   = find_bits(p.rxChanMask);
    if isempty(rxList), rxList = 1:numRx; end
    rxRef    = rxList(1);

    if isTDM
        h_rx = mean(Htxrx, 1);   % average across TXs, size [1 x numRx]
    else
        h_rx = Htxrx(1,:);       % only RX info available
    end
    % make reference channel real-positive; equalize others to it
    refVal  = h_rx(rxRef);
    rxW     = conj(h_rx ./ refVal);
    rxW(~isfinite(rxW)) = 0;

    % ---- TX weights (only if TDM) ----
    if isTDM
        % after RX equalization, compute per-TX phase wrt reference TX
        txRefIdx = 1;                % choose the first active TX as reference
        h_tx = zeros(numTxAct,1);
        for tIdx = 1:numTxAct
            h_tx(tIdx) = sum( Htxrx(tIdx,:) .* rxW );  % combine across RX after RX EQ
        end
        txRefVal = h_tx(txRefIdx);
        txW = conj(h_tx ./ txRefVal);
        txW(~isfinite(txW)) = 0;
    else
        txW = [];
    end

    % ---- range bias (measured - known) ----
    R_meas = k0 * dR;
    if isempty(o.cornerRange_m)
        rangeBias = NaN;
    else
        rangeBias = R_meas - o.cornerRange_m;
    end

    % ---- pack results ----
    cal = struct();
    cal.rxWeights    = rxW(:);
    cal.txWeights    = txW(:);
    cal.rangeBias_m  = rangeBias;
    cal.rangeBin     = k0;
    cal.dR_m         = dR;
    cal.kdopp        = kdopp;
    cal.meta = struct( ...
        'isTDM', isTDM, ...
        'txActive', txActive(:), ...
        'patternLen', patternLen, ...
        'txOfChirpPattern', txOfChirp(1:patternLen), ...
        'Fs_Hz', Fs, 'Slope_Hz_per_s', S, ...
        'NfftR', o.NfftR, 'NfftD', o.NfftD, ...
        'framesUsed', useF, ...
        'rxRef', rxRef, ...
        'R_meas_m', R_meas );

    % ---- quick sanity plots (optional, robust to shapes) ----
    if o.plot
        figure('Name','Corner calibration','Color','w');

        % Subplot 1: Range spectrum around search window with vertical line at R_meas
        subplot(2,2,1);
        rAx = (0:o.NfftR-1)*dR;             % meters
        x = rAx(rSel);
        y = P_r(:);
        if numel(y) ~= numel(x)
            n = min(numel(y), numel(x));
            x = x(1:n); y = y(1:n);
        end
        plot(x, y); grid on; hold on;
        xline(R_meas,'r--','HandleVisibility','off');
        title('Range spectrum (sum over RX/chirps)');
        xlabel('Range (m)'); ylabel('Power');

        % Subplot 2: Per-RX phase at corner (pre-EQ)
        subplot(2,2,2);
        stem(1:numRx, angle(h_rx)*180/pi, 'filled'); grid on;
        xlim([0.5 numRx+0.5]); xticks(1:numRx);
        title('Per-RX phase at corner (pre-EQ)'); ylabel('deg'); xlabel('RX index');

        % Subplot 3: Per-RX magnitude at corner (pre-EQ)
        subplot(2,2,3);
        mm = abs(h_rx); mm = mm / max(eps, max(mm));
        stem(1:numRx, mm, 'filled'); grid on;
        xlim([0.5 numRx+0.5]); xticks(1:numRx);
        title('Per-RX magnitude at corner (pre-EQ)'); ylabel('norm mag'); xlabel('RX index');

        % Subplot 4: RX equalization weights (phase)
        subplot(2,2,4);
        stem(1:numRx, angle(rxW)*180/pi, 'filled'); grid on;
        xlim([0.5 numRx+0.5]); xticks(1:numRx);
        title('RX equalization weights (phase)'); ylabel('deg'); xlabel('RX index');
    end
end

% ===== helpers =====

function [txOfChirp, patternLen, txActive, isTDM] = tx_mapping_from_params(p)
% Build per-chirp TX assignment within one "pattern"
% Returns:
%   txOfChirp : [1 x patternLen] values in {0,1,2,3} (0=unknown/multi), 1..3 for TX1..TX3
%   patternLen: number of distinct chirp IDs in a pattern
%   txActive  : list of active TX indices (1..3)
%   isTDM     : true if every chirp enables exactly one TX (power-of-two mask)

    if ~isfield(p,'chirpSegments') || isempty(p.chirpSegments)
        % Fall back: single-chirp pattern with unknown TX
        txOfChirp = 0; patternLen = 1; txActive = find_bits(p.txChanMask);
        isTDM = (numel(txActive)==1);
        return;
    end
    seg = p.chirpSegments; % rows: [start end txMask]
    chirpList = [];
    maskList  = [];
    for i=1:size(seg,1)
        cId = seg(i,1):seg(i,2);
        chirpList = [chirpList, cId];
        maskList  = [maskList, repmat(seg(i,3), 1, numel(cId))];
    end
    [chirpList, ord] = sort(chirpList);
    maskList = maskList(ord);
    patternLen = numel(chirpList);

    txOfChirp = zeros(1, patternLen);
    txActive = [];
    for k=1:patternLen
        m = maskList(k);
        if is_power_of_two(m) && m~=0
            txOfChirp(k) = log2(double(m)) + 1;  % 1..3
        else
            txOfChirp(k) = 0;                     % multi-TX or unknown
        end
        txActive = union(txActive, find_bits(m));
    end
    isTDM = all(txOfChirp~=0) && ~isempty(txActive);
end

function tf = is_power_of_two(x)
    x = uint32(x);
    tf = (x~=0) && bitand(x, x-1)==0;
end

function idx = find_bits(mask)
% Return 1-based indices set in mask (LSB=1)
    idx = [];
    for k=1:32
        if bitand(uint32(mask), bitshift(uint32(1),k-1))
            idx(end+1) = k; %#ok<AGROW>
        end
    end
end
