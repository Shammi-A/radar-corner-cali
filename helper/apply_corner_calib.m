function cubeCal = apply_corner_calib(cube, p, cal)
% APPLY_CORNER_CALIB  Apply RX/TX equalization weights to raw cube
% cube   : [numRx x Ns x Nc_pf x Nf]
% p      : RadarParams (for chirp mapping)
% cal    : from corner_calibrate (rxWeights, txWeights, meta)
%
% Returns calibrated cube of the same size.

    cubeCal = cube;

    % RX weights (always)
    if isfield(cal,'rxWeights') && ~isempty(cal.rxWeights)
        Wrx = reshape(cal.rxWeights(:), [], 1, 1, 1);  % [numRx x 1 x 1 x 1]
        cubeCal = bsxfun(@times, cubeCal, Wrx);
    end

    % TX weights (only if TDM)
    if isfield(cal,'txWeights') && ~isempty(cal.txWeights) && isfield(cal,'meta') && cal.meta.isTDM
        [txOfChirp, patternLen, txActive] = tx_mapping_from_params_local(p);
        loops = max(1, p.chirpsPerLoop);
        Nc_pf = size(cube,3);

        % Expand per-chirp TX indices across the frame
        txCh = zeros(1, Nc_pf);
        % pattern repeats 'loops' times
        for l = 1:loops
            s = (l-1)*patternLen + 1;
            e = min(Nc_pf, l*patternLen);
            L = e - s + 1;
            txCh(s:e) = txOfChirp(1:L);
        end

        % Apply per-chirp TX weight (multiply across RX for those chirps)
        for tIdx = 1:numel(txActive)
            tx = txActive(tIdx);
            w  = cal.txWeights(tIdx);
            chirpIdx = find(txCh == tx);
            if isempty(chirpIdx), continue; end
            cubeCal(:,:,chirpIdx,:) = cubeCal(:,:,chirpIdx,:) * w;
        end
    end
end

% local copy to avoid external dependency
function [txOfChirp, patternLen, txActive] = tx_mapping_from_params_local(p)
    if ~isfield(p,'chirpSegments') || isempty(p.chirpSegments)
        txOfChirp = 0; patternLen = 1; txActive = find_bits(p.txChanMask); return;
    end
    seg = p.chirpSegments;
    chirpList = []; maskList = [];
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
            txOfChirp(k) = 0;                     % multi-TX/unknown → cannot TX-cal
        end
        txActive = union(txActive, find_bits(m));
    end
end

function tf = is_power_of_two(x)
    x = uint32(x);
    tf = (x~=0) && bitand(x, x-1)==0;
end

function idx = find_bits(mask)
    idx = [];
    for k=1:3
        if bitand(uint32(mask), bitshift(uint32(1),k-1))
            idx(end+1) = k; %#ok<AGROW>
        end
    end
end
