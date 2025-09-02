function [RadarParams, Raw] = mmws_parse_log(logFile, opts)
% MMWS_PARSE_LOG  Robust TI mmWave Studio log parser (AWR2243 + DCA1000)
% - Exact-key parsing for API lines (no substring collisions)
% - TX/RX masks in any order (incl. hex)
% - Slope from device codeword using band-appropriate LSB (77/60 GHz)
% - Profile times auto-units (µs or 10-ns ticks)
% - FrameConfig (preferred) or AdvancedFrameConfig; 5-ns LSB for large period values
% - ChirpConfig → chirp segments + chirps per frame (loops × segments)
% - Doppler metrics: bins, resolution, v_max, frame rate

    if nargin < 2
        opts = struct();
    end
    if ~isfield(opts,'profileId')
        opts.profileId = [];
    end

    assert(exist(logFile,'file')==2, 'File not found: %s', logFile);

    % --- read file ---
    L = read_lines(logFile);
    whole = strjoin(L, newline);

    % --- helpers to fetch args (EXACT key match only) ---
    get1   = @(k) get_args_one_exact(L, k);
    getAll = @(k) get_args_all_exact(L, k);

    % --- capture raw (last occurrence) ---
    Raw = struct();
    Raw.ChannelConfig        = first_nonempty({get1('ChannelConfig'),        get1('ar1.ChannelConfig')});
    Raw.ProfileConfig        = first_nonempty({get1('ProfileConfig'),        get1('ar1.ProfileConfig')});
    Raw.ChirpConfig          = first_nonempty({get1('ChirpConfig'),          get1('ar1.ChirpConfig')});
    Raw.FrameConfig          = first_nonempty({get1('FrameConfig'),          get1('ar1.FrameConfig')});
    Raw.AdvancedFrameConfig  = first_nonempty({get1('AdvancedFrameConfig'),  get1('ar1.AdvanceFrameConfig')});

    % also ALL occurrences (for multi-profile / multi-chirp)
    Raw.ProfileConfig_ALL    = [getAll('ProfileConfig'), getAll('ar1.ProfileConfig')];
    Raw.ChirpConfig_ALL      = [getAll('ChirpConfig'),   getAll('ar1.ChirpConfig')];

    % optional extras (unused but parsed for completeness)
    Raw.AdcOutConfig         = first_nonempty({get1('AdcOutConfig'),         get1('ar1.ADCOutCfg')});
    Raw.DataFmtConfig        = first_nonempty({get1('DataFmtConfig'),        get1('ar1.DataFmtConfig')});

    % --- HSI/TSW hints (if present) ---
    [hsiHz, tswSamp] = parse_tsw(L);

    % --- TX/RX masks ---
    [txMask, rxMask] = interpret_masks(Raw.ChannelConfig);
    p.txChanMask = txMask;
    p.rxChanMask = rxMask;
    p.numTx = bitcount_u32(txMask);
    p.numRx = bitcount_u32(rxMask);

    % --- choose profile (if multiple) ---
    pc = choose_profile(Raw, opts.profileId);

    % TI ProfileConfig expected order:
    % [profileId, startFreqConst, idle, adcStart, rampEnd, txPow, txPhase,
    %  slopeCodeword, txStart, numAdcSamples, digOutRate_ksps, ...]
    startFreqConst = geti(pc,2);
    idle           = geti(pc,3);
    adcStart       = geti(pc,4);
    rampEnd        = geti(pc,5);
    slopeCodeword  = geti(pc,8);
    txStart        = geti(pc,9);
    numAdc         = geti(pc,10);
    digKsps        = geti(pc,11);  % ksps

    % --- slope from device codeword (band-aware LSB) ---
    % Use startFreqConst to decide band; ~53.644 Hz/LSB synth codeword.
    f0_Hz = double(startFreqConst) * 53.644;
    if isfinite(f0_Hz) && f0_Hz > 70e9
        lsb_kHz_per_us = 48.279;                   % 77 GHz family (AWR2243)
    else
        lsb_kHz_per_us = 36.210;                   % 60 GHz family
    end
    p.slope_MHz_per_us = slopeCodeword * (lsb_kHz_per_us / 1000.0);  % MHz/us
    p.Slope_Hz_per_s   = p.slope_MHz_per_us * 1e12;                  % Hz/s

    % --- profile time units (µs vs 10-ns ticks) ---
    u_prof = detect_profile_unit([idle adcStart rampEnd txStart]);   % 'us' or 'tick10ns'
    p.idleTime_s     = conv_time(idle,     u_prof);
    p.adcStartTime_s = conv_time(adcStart, u_prof);
    p.rampEndTime_s  = conv_time(rampEnd,  u_prof);
    p.txStartTime_s  = conv_time(txStart,  u_prof);
    p.T_chirp_s      = p.idleTime_s + p.rampEndTime_s;

    % --- sampling / ADC ---
    p.numADCSamples  = numAdc;
    p.Fs_Hz          = digKsps * 1e3;

    % --- bandwidths & range resolution ---
    p.B_total_Hz     = p.Slope_Hz_per_s * nanzero(p.rampEndTime_s);       % full ramp BW (GUI)
    if all(~isnan([p.numADCSamples, p.Fs_Hz])) && p.Fs_Hz > 0
        T_adc = p.numADCSamples / p.Fs_Hz;
        p.B_sampled_Hz = p.Slope_Hz_per_s * T_adc;                        % BW seen by range FFT
        p.rangeRes_m   = 3e8 / (2 * p.B_sampled_Hz);
    else
        p.B_sampled_Hz = NaN;
        p.rangeRes_m   = NaN;
    end

    % --- frame: prefer basic FrameConfig; else Advanced ---
    if ~isempty(Raw.FrameConfig)
        [numFrames, loops, framePeriod_s] = parse_frame(Raw.FrameConfig);
    elseif ~isempty(Raw.AdvancedFrameConfig)
        [numFrames, loops, framePeriod_s] = parse_adv_frame(Raw.AdvancedFrameConfig);
    else
        numFrames = NaN; loops = NaN; framePeriod_s = NaN;
    end
    % Fallback to labeled text if needed
    if isnan(loops) || isnan(numFrames) || isnan(framePeriod_s)
        [loops2, frames2, per2] = parse_frame_from_labels(whole);
        if isnan(loops);         loops = loops2;       end
        if isnan(numFrames);     numFrames = frames2;  end
        if isnan(framePeriod_s); framePeriod_s = per2; end
    end
    if isnan(loops); loops = 1; end
    if isnan(numFrames); numFrames = 1; end

    p.numFrames     = numFrames;
    p.chirpsPerLoop = loops;                 % "loops" = chirps per TX per frame
    p.framePeriod_s = framePeriod_s;
    p.frameRate_Hz  = 1./p.framePeriod_s;

    % --- chirp segments & chirps per frame ---
    [segments, nChirpsOnce] = summarize_chirps_robust(Raw.ChirpConfig_ALL, pc);
    p.chirpSegments  = segments;                             % [startCh, endCh, txMask]
    p.chirpsPerFrame = max(1, nChirpsOnce) * max(1, round(loops));

    % --- Doppler metrics (TDM-aware & conservative) ---
    lambda = compute_lambda(f0_Hz);
    nChirpsPerPattern = 0;
    if ~isempty(p.chirpSegments)
        nChirpsPerPattern = sum(max(0, p.chirpSegments(:,2) - p.chirpSegments(:,1) + 1));
    end
    if ~isempty(p.chirpSegments)
        txBits = arrayfun(@(m) sum(dec2bin(uint32(m))=='1'), p.chirpSegments(:,3));
        nTxPerPattern = max(1, max(txBits));
    else
        nTxPerPattern = max(1, p.numTx);
    end
    if nChirpsPerPattern >= nTxPerPattern && nTxPerPattern > 0
        Ndopp = max(1, floor(nChirpsPerPattern / nTxPerPattern) * round(p.chirpsPerLoop));
    else
        Ndopp = max(1, round(p.chirpsPerLoop));  % safe fallback
    end
    p.dopplerBins  = Ndopp;
    p.vRes_mps     = lambda / (2 * p.T_chirp_s * Ndopp);
    p.vMax_mps     = lambda / (4 * p.T_chirp_s);             % ±vMax unambiguous

    % --- extras ---
    p.fc_const_raw        = startFreqConst;
    p.hsiClk_Hz_fromLine  = hsiHz;
    p.tsw_sampRate_val    = tswSamp;
    p.lambda_m            = lambda;
    p.units_profile_times = u_prof;
    p.notes = 'Exact-key parse; slope from codeword LSB; FrameConfig(5ns); robust TX mask; Doppler metrics.';

    RadarParams = p;
end

% ================= helpers =================

function L = read_lines(f)
    raw = fileread(f);
    raw = strrep(raw, sprintf('\r\n'), sprintf('\n'));
    raw = strrep(raw, sprintf('\r'),   sprintf('\n'));
    L = regexp(raw, '\n', 'split')';
    L = L(~cellfun(@isempty,L));
end

function v = get_args_one_exact(L, key)
% Return args for the LAST line that is exactly "API:key" or "ar1.key("
    v = [];
    pat_api = ['API:\s*' regexptranslate('escape', key) '\s*,'];
    pat_lua = ['ar1\.'   regexptranslate('escape', key) '\s*\('];
    hits = false(numel(L),1);
    for i=1:numel(L)
        s = L{i};
        if ~isempty(regexp(s, pat_api, 'once')) || ~isempty(regexp(s, pat_lua, 'once'))
            hits(i) = true;
        end
    end
    idx = find(hits);
    if isempty(idx), return; end
    S = L{idx(end)};
    m = regexp(S, [regexptranslate('escape',key) '\s*\(([^)]*)\)'], 'tokens', 'once');  % ar1.key(...)
    if ~isempty(m)
        v = parse_args(m{1}); return;
    end
    m = regexp(S, ['API:\s*' regexptranslate('escape',key) '\s*,(.*)$'], 'tokens', 'once'); % API:key, ...
    if ~isempty(m)
        v = parse_args(m{1}); return;
    end
end

function A = get_args_all_exact(L, key)
% Return args for ALL lines that are exactly "API:key" or "ar1.key("
    A = {};
    pat_api = ['API:\s*' regexptranslate('escape', key) '\s*,'];
    pat_lua = ['ar1\.'   regexptranslate('escape', key) '\s*\('];
    for i=1:numel(L)
        s = L{i};
        if isempty(regexp(s, pat_api, 'once')) && isempty(regexp(s, pat_lua, 'once'))
            continue;
        end
        m = regexp(s, [regexptranslate('escape',key) '\s*\(([^)]*)\)'], 'tokens', 'once');
        if isempty(m)
            m = regexp(s, ['API:\s*' regexptranslate('escape',key) '\s*,(.*)$'], 'tokens', 'once');
        end
        if ~isempty(m)
            A{end+1} = parse_args(m{1}); %#ok<AGROW>
        end
    end
end

function nums = parse_args(argstr)
    toks = regexp(argstr, ',', 'split');
    toks = strtrim(toks);
    toks = toks(~cellfun(@isempty,toks));
    nums = zeros(1,numel(toks));
    for i = 1:numel(toks)
        t = regexprep(toks{i}, '[;)]$', '');
        if ~isempty(regexp(t, '^0x[0-9A-Fa-f]+$', 'once'))
            nums(i) = double(hex2dec(t(3:end)));
        elseif ~isempty(regexp(t, '^[\+\-]?\d+(\.\d+)?([eE][\+\-]?\d+)?$', 'once'))
            nums(i) = str2double(t);
        else
            nums(i) = NaN; % ignore strings etc.
        end
    end
    while ~isempty(nums) && isnan(nums(end))
        nums(end) = [];
    end
end

function out = first_nonempty(list)
    out = [];
    for i = 1:numel(list)
        if ~isempty(list{i})
            out = list{i};
            return;
        end
    end
end

function [hsiHz, tswSamp] = parse_tsw(L)
    hsiHz = NaN; tswSamp = NaN;
    idx = find(contains(L, 'Sampling rate'));
    if isempty(idx)
        return;
    end
    S = L{idx(end)};
    nums = regexp(S, '\d+', 'match');
    if numel(nums) >= 2
        hsiHz  = str2double(nums{1});
        tswSamp= str2double(nums{2});
    end
end

function [txMask, rxMask] = interpret_masks(v)
    txMask = uint32(0);
    rxMask = uint32(0);
    if isempty(v) || numel(v) < 2
        return;
    end
    a = sanitize_mask(v(1));
    b = sanitize_mask(v(2));
    if bitcount_u32(a) <= 3 && bitcount_u32(b) <= 4
        txMask = a; rxMask = b;
    elseif bitcount_u32(b) <= 3 && bitcount_u32(a) <= 4
        txMask = b; rxMask = a;
    else
        if a <= b
            txMask = a; rxMask = b;
        else
            txMask = b; rxMask = a;
        end
    end
end

function u = sanitize_mask(x)
    if isnan(x)
        u = uint32(0); return;
    end
    x = floor(abs(x));
    u = uint32(mod(x, 2^32));
end

function n = bitcount_u32(x)
    n = sum(dec2bin(double(x))=='1');
end

function y = geti(v,i)
    if numel(v) >= i
        y = v(i);
    else
        y = NaN;
    end
end

function unit = detect_profile_unit(vals)
    vals = vals(~isnan(vals));
    if isempty(vals)
        unit = 'us';
        return;
    end
    % 10-ns ticks often yield values like 10000 for 100 µs → median >= 2000
    if median(vals) >= 2000
        unit = 'tick10ns';
    else
        unit = 'us';
    end
end

function t = conv_time(x, unit)
    if isnan(x)
        t = NaN; return;
    end
    switch unit
        case 'us'
            t = x * 1e-6;
        case 'tick10ns'
            t = x * 10e-9;
        otherwise
            t = NaN;
    end
end

function z = nanzero(x)
    if isnan(x) || x == 0
        z = eps;
    else
        z = x;
    end
end

function [numFrames, loops, framePeriod_s] = parse_frame(fc)
% ar1.FrameConfig(startChirp, endChirp, numFrames, numLoops, framePeriod, ...)
% NOTE: Studio log order = frames, loops (not SDK CLI order)
% Period typically in 5-ns ticks when large (e.g., 8,000,000 -> 40 ms)
    numFrames = NaN; loops = NaN; framePeriod_s = NaN;
    if isempty(fc) || numel(fc) < 5
        return;
    end
    numFrames = fc(3);
    loops     = fc(4);
    per       = fc(5);
    if per > 1e5
        framePeriod_s = per * 5e-9;      % 5-ns LSB
    else
        framePeriod_s = per * 1e-3;      % ms (if Studio ever prints ms)
    end
end

function [numFrames, loops, framePeriod_s] = parse_adv_frame(af)
% Best-effort for AdvancedFrameConfig (period in ms or ticks somewhere)
    numFrames = NaN; loops = NaN; framePeriod_s = NaN;
    if isempty(af)
        return;
    end
    ms  = af(af >= 1   & af <= 2000);
    tks = af(af >= 5e5 & af <= 2e8);
    if ~isempty(ms)
        framePeriod_s = ms(1) * 1e-3;
    elseif ~isempty(tks)
        framePeriod_s = tks(1) * 5e-9;   % most Studio builds use 5-ns ticks here too
    end
    ints = af(af > 0 & af < 1e7 & abs(af - round(af)) < 1e-6);
    if ~isempty(ints)
        numFrames = ints(1);
    end
end

function [loops, frames, per_s] = parse_frame_from_labels(text)
% Fallback: scan labeled lines anywhere in the log
    loops  = NaN; frames = NaN; per_s = NaN;
    m = regexp(text, 'No of Chirp Loops[^0-9]*(\d+)', 'tokens', 'once');
    if ~isempty(m); loops = str2double(m{1}); end
    m = regexp(text, 'No of Frames[^0-9]*(\d+)', 'tokens', 'once');
    if ~isempty(m); frames = str2double(m{1}); end
    m = regexp(text, 'Periodicity\s*\(ms\)[^0-9]*(\d+(\.\d+)?)', 'tokens', 'once');
    if ~isempty(m); per_s = str2double(m{1}) * 1e-3; end
end

function [segments, nChirpsOnce] = summarize_chirps_robust(ChList, pc)
% Build [startCh, endCh, txMask] rows; count chirps per once-through for the chosen profile
    segments = [];
    nChirpsOnce = 0;
    if isempty(ChList)
        return;
    end
    profId = pc(1);
    for k = 1:numel(ChList)
        v = ChList{k};
        if numel(v) < 3
            continue;
        end
        if v(1) ~= profId
            continue;
        end
        startCh = max(0, round(v(2)));
        endCh   = max(startCh, round(v(3)));
        txMask  = infer_tx_mask_from_tail(v);
        segments = [segments; startCh, endCh, txMask]; %#ok<AGROW>
        nChirpsOnce = nChirpsOnce + (endCh - startCh + 1);
    end
end

function txMask = infer_tx_mask_from_tail(v)
% Tolerate two formats:
%  A) single TX bitmask near the end (e.g., ..., 7, 0)
%  B) trailing three booleans (TX0, TX1, TX2)
    txMask = 0;
    if numel(v) <= 7
        return;
    end
    tail = v(8:end);

    % Case A: last integer in [1..7] ⇒ mask
    ints = tail(~isnan(tail) & abs(tail - round(tail)) < 1e-6);
    idxMask = find(ints >= 1 & ints <= 7, 1, 'last');
    if ~isempty(idxMask)
        txMask = uint32(bitand(ints(idxMask), 7));
        return;
    end

    % Case B: last three boolean-ish entries ⇒ [TX0 TX1 TX2]
    boolIdx = find(~isnan(tail) & (tail==0 | tail==1));
    if numel(boolIdx) >= 3
        b = tail(boolIdx(end-2:end));
        for k = 1:3
            if b(k) == 1
                txMask = bitor(txMask, bitshift(1, k-1));
            end
        end
        return;
    end

    % Fallback: last 3 entries nonzero → set bits
    if numel(tail) >= 3
        b = double(tail(end-2:end) ~= 0);
        for k = 1:3
            if b(k)
                txMask = bitor(txMask, bitshift(1, k-1));
            end
        end
    end
end

function pc = choose_profile(Raw, preferId)
% Choose a ProfileConfig by preferred ID; else prefer id==0; else last one.
    P = Raw.ProfileConfig_ALL;
    if isempty(P)
        pc = Raw.ProfileConfig;
        return;
    end
    if ~isempty(preferId)
        idx = [];
        for k = 1:numel(P)
            if ~isempty(P{k}) && P{k}(1) == preferId
                idx(end+1) = k; %#ok<AGROW>
            end
        end
        if ~isempty(idx)
            pc = P{idx(end)};
            return;
        end
    end
    idx0 = [];
    for k = 1:numel(P)
        if ~isempty(P{k}) && P{k}(1) == 0
            idx0(end+1) = k; %#ok<AGROW>
        end
    end
    if ~isempty(idx0)
        pc = P{idx0(end)};
    else
        pc = P{end};
    end
end

function lambda = compute_lambda(f0_Hz)
% Use start frequency if reasonable; else default to 77 GHz
    if isfinite(f0_Hz) && f0_Hz > 50e9 && f0_Hz < 90e9
        lambda = 3e8 / f0_Hz;
    else
        lambda = 3e8 / 77e9;
    end
end
