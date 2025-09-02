function stats = corner_calib_check(cube, p, cal, varargin)
% CORNER_CALIB_CHECK  One-click validation of corner-reflector calibration.
% Inputs
%   cube : [numRx x Ns x Nc_pf x Nf] raw cube (from dca1000_read_bin)
%   p    : RadarParams (from mmws_parse_log)
%   cal  : struct from corner_calibrate (rxWeights, txWeights, rangeBin, dR_m, ...)
%
% Name-Value options (all optional):
%   'useFrames'     : how many frames to average (default: min(4, Nf))
%   'cornerRange_m' : known distance (overrides cal.rangeBin if given)
%   'NfftR'         : range FFT size (default: nextpow2(2*Ns))
%   'plot'          : true/false (default: true)
%
% Output 'stats' fields:
%   .R_est_m, .R_known_m, .rangeBias_m
%   .SNR_before_dB, .SNR_after_dB, .SNR_gain_dB
%   .phaseSpread_before_deg, .phaseSpread_after_deg, .phaseSpread_gain_deg
%   .magImbalance_before_dB, .magImbalance_after_dB, .magImbalance_gain_dB

    ip = inputParser;
    ip.addParameter('useFrames', []);
    ip.addParameter('cornerRange_m', []);
    ip.addParameter('NfftR', []);
    ip.addParameter('plot', true);
    ip.parse(varargin{:});
    o = ip.Results;

    [numRx, Ns, Nc_pf, Nf] = size(cube);
    useF = min(max(1, ~isempty(o.useFrames) * o.useFrames + isempty(o.useFrames) * min(4,Nf)), Nf);

    if isempty(o.NfftR), o.NfftR = 2^nextpow2(max(512, 2*Ns)); end

    % -------- Range FFTs (before & after) --------
    % Average 'useF' frames coherently across frames, incoherently over chirps
    wR = hann(Ns,'periodic').'; wR = wR / max(eps, norm(wR));

    % Before calibration
    X  = reshape(cube(:,:,:,1:useF), numRx, Ns, Nc_pf*useF);
    Xw = bsxfun(@times, X, reshape(wR,1,Ns,1));
    Xr = fft(Xw, o.NfftR, 2);                         % [numRx x NfftR x (Nc_pf*useF)]
    P_before = squeeze(sum(abs(Xr).^2, [1 3]));       % [NfftR x 1], summed over RX & shots

    % After calibration
    cubeCal = apply_corner_calib(cube(:,:,:,1:useF), p, cal);
    Xc  = reshape(cubeCal, numRx, Ns, Nc_pf*useF);
    Xcw = bsxfun(@times, Xc, reshape(wR,1,Ns,1));
    Xcr = fft(Xcw, o.NfftR, 2);
    P_after = squeeze(sum(abs(Xcr).^2, [1 3]));

    % -------- Pick corner bin (from known range or from cal / peak) --------
    dR = 3e8 * (p.Fs_Hz/o.NfftR) / (2 * p.Slope_Hz_per_s);
    if ~isempty(o.cornerRange_m)
        kCorner = max(1, min(floor(o.NfftR/2), round(o.cornerRange_m/dR)));
        R_known = o.cornerRange_m;
    elseif isfield(cal,'rangeBin') && ~isempty(cal.rangeBin) && isfinite(cal.rangeBin)
        kCorner = min(cal.rangeBin, floor(o.NfftR/2));
        R_known = isfield(cal,'rangeBias_m') * (cal.rangeBin*cal.dR_m - cal.rangeBias_m) + ~isfield(cal,'rangeBias_m')*NaN;
    else
        % fall back: use the global peak in the positive half (before-cal)
        [~, kCorner] = max(P_before(1:floor(o.NfftR/2)));
        R_known = NaN;
    end
    R_est = kCorner * dR;

    % -------- SNR around the corner peak --------
    % Exclude ±2 bins around the peak for noise estimate
    pos = 1:floor(o.NfftR/2);
    excl = max(1,kCorner-2):min(floor(o.NfftR/2),kCorner+2);
    noiseIdx = setdiff(pos, excl);

    SNR_b = 10*log10( (P_before(kCorner)+eps) / (median(P_before(noiseIdx))+eps) );
    SNR_a = 10*log10( (P_after(kCorner)+eps)  / (median(P_after(noiseIdx))+eps)  );

    % -------- Per-RX phase & magnitude at the corner (before & after) --------
    % Slice corner bin across shots, then average to get one complex h per RX
    Xk_b = reshape(Xr(:,kCorner,:), numRx, Nc_pf*useF);   % before
    Xk_a = reshape(Xcr(:,kCorner,:), numRx, Nc_pf*useF);  % after
    h_b  = mean(Xk_b, 2);                                  % [numRx x 1]
    h_a  = mean(Xk_a, 2);

    phSpread_b = rad2deg(ang_spread(h_b));
    phSpread_a = rad2deg(ang_spread(h_a));

    magImb_b_dB = 20*log10( max(abs(h_b))/max(eps, min(abs(h_b))) );
    magImb_a_dB = 20*log10( max(abs(h_a))/max(eps, min(abs(h_a))) );

    % -------- Report --------
    stats = struct();
    stats.R_est_m  = R_est;
    stats.R_known_m = R_known;
    stats.rangeBias_m = R_est - R_known;  % NaN if R_known is NaN

    stats.SNR_before_dB = SNR_b;
    stats.SNR_after_dB  = SNR_a;
    stats.SNR_gain_dB   = SNR_a - SNR_b;

    stats.phaseSpread_before_deg = phSpread_b;
    stats.phaseSpread_after_deg  = phSpread_a;
    stats.phaseSpread_gain_deg   = phSpread_b - phSpread_a;

    stats.magImbalance_before_dB = magImb_b_dB;
    stats.magImbalance_after_dB  = magImb_a_dB;
    stats.magImbalance_gain_dB   = magImb_b_dB - magImb_a_dB;

    fprintf('\n=== Corner calibration check ===\n');
    fprintf('Range (est) = %.3f m', R_est);
    if ~isnan(R_known), fprintf('   | known = %.3f m   | bias = %.3f m', R_known, stats.rangeBias_m); end
    fprintf('\n');
    fprintf('SNR: before = %6.2f dB   after = %6.2f dB   gain = %+6.2f dB\n', SNR_b, SNR_a, stats.SNR_gain_dB);
    fprintf('RX phase spread: before = %6.2f°  after = %6.2f°  improvement = %+6.2f°\n', phSpread_b, phSpread_a, stats.phaseSpread_gain_deg);
    fprintf('RX mag imbalance: before = %5.2f dB  after = %5.2f dB  improvement = %+5.2f dB\n', magImb_b_dB, magImb_a_dB, stats.magImbalance_gain_dB);

    % -------- Plots --------
    if o.plot
        N2 = floor(o.NfftR/2);
        ra = (0:o.NfftR-1)*dR;
        figure('Name','Corner calib check','Color','w');
        subplot(1,2,1);
        plot(ra(1:N2), 10*log10(P_before(1:N2)+eps), 'DisplayName','before'); hold on;
        plot(ra(1:N2), 10*log10(P_after(1:N2)+eps),  'DisplayName','after');
        xline(R_est,'r--','HandleVisibility','off'); grid on; legend location best;
        xlabel('Range (m)'); ylabel('Power (dB)'); title('Range spectrum @ sum(RX,shots)');

        subplot(1,2,2);
        stem(1:numRx, angle(h_b)*180/pi, 'o','DisplayName','before'); hold on;
        stem(1:numRx, angle(h_a)*180/pi, 'filled','DisplayName','after'); grid on; legend;
        xlim([0.5 numRx+0.5]); xticks(1:numRx);
        xlabel('RX index'); ylabel('deg'); title('Per-RX phase @ corner');
    end
end

% ----- helpers -----
function s = ang_spread(h)
% Circular standard deviation of phases of complex vector h
    a = angle(h(:));
    c = mean(exp(1j*a));
    s = sqrt(2*(1-abs(c)));   % small-angle approx for reporting
end
