function [RadarParams, Raw] = mmws_parse_log(logFile, opts)
% Robust TI mmWave Studio log parser (AWR2243 + DCA1000)
% - Args-only extraction for API lines (Lua & "API:" styles)
% - TX/RX masks in any order (incl. hex)
% - Slope unit normalization: kHz/us | 0.1 MHz/us | MHz/us
% - Profile times: auto-units (us | 10-ns ticks)
% - FrameConfig: also parsed from *text labels* if ar1.FrameConfig(...) absent
% - ChirpConfig: robust TX mask from trailing boolean flags (TX0/TX1/TX2)

    if nargin < 2, opts = struct(); end
    if ~isfield(opts,'profileId'), opts.profileId = []; end

    assert(exist(logFile,'file')==2, 'File not found: %s', logFile);
    L = read_lines(logFile);                 % cellstr lines
    whole = strjoin(L, newline);            % full text for label fallback

    % --- shorthands ---
    get1   = @(k) get_args_one(L, k);
    getAll = @(k) get_args_all(L, k);

    % --- raw captures (last occurrence) ---
    Raw = struct();
    Raw.ChannelConfig        = first_nonempty({get1('ChannelConfig'),        get1('ar1.ChannelConfig')});
    Raw.ProfileConfig        = first_nonempty({get1('ProfileConfig'),        get1('ar1.ProfileConfig')});
    Raw.ChirpConfig          = first_nonempty({get1('ChirpConfig'),          get1('ar1.ChirpConfig')});
    Raw.FrameConfig          = first_nonempty({get1('FrameConfig'),          get1('ar1.FrameConfig')});
    Raw.AdvancedFrameConfig  = first_nonempty({get1('AdvancedFrameConfig'),  get1('ar1.AdvanceFrameConfig')});

    % also collect ALL occurrences for multi-profile/multi-chirp
    Raw.ProfileConfig_ALL    = [getAll('ProfileConfig'),      getAll('ar1.ProfileConfig')];
    Raw.ChirpConfig_ALL      = [getAll('ChirpConfig'),        getAll('ar1.ChirpConfig')];

    % optional extras (unused here but left for completeness)
    Raw.AdcOutConfig         = first_nonempty({get1('AdcOutConfig'),         get1('ar1.ADCOutCfg')});
    Raw.DataFmtConfig        = first_nonempty({get1('DataFmtConfig'),        get1('ar1.DataFmtConfig')});

    % --- HSI/TSW hint (if present) ---
    [hsiHz, tswSamp] = parse_tsw(L);

    % --- Channel masks ---
    [txMask, rxMask] = interpret_masks(Raw.ChannelConfig);
    p.txChanMask = txMask;             p.rxChanMask = rxMask;
    p.numTx = bitcount_u32(txMask);    p.numRx = bitcount_u32(rxMask);

    % --- Choose Profile (if multiple) ---
    pc = choose_profile(Raw, opts.profileId);

    % expected TI order:
    % [profileId, startFreqConst, idleTime, adcStartTime, rampEndTime,
    %  txOutPow, txPhase, freqSlopeConst, txStartTime, numAdc, digOutRate(ksps), ...]
    startFreqConst = geti(pc,2);
    idle           = geti(pc,3);
    adcStart       = geti(pc,4);
    rampEnd        = geti(pc,5);
    rawSlope       = geti(pc,8);
    txStart        = geti(pc,9);
    numAdc         = geti(pc,10);
    digKsps        = geti(pc,11);

    % --- Slope normalization ---
    if     rawSlope > 1000, p.slope_MHz_per_us = rawSlope/1000;    % kHz/us
    elseif rawSlope > 100,  p.slope_MHz_per_us = rawSlope/10;      % 0.1 MHz/us
    else                    p.slope_MHz_per_us = rawSlope;         % MHz/us
    end
    p.Slope_Hz_per_s = p.slope_MHz_per_us * 1e12;

    % --- Profile time units ---
    u_prof = detect_profile_unit([idle adcStart rampEnd txStart]);  % 'us'|'tick10ns'
    p.idleTime_s     = conv_time(idle,     u_prof);
    p.adcStartTime_s = conv_time(adcStart, u_prof);
    p.rampEndTime_s  = conv_time(rampEnd,  u_prof);
    p.txStartTime_s  = conv_time(txStart,  u_prof);
    p.T_chirp_s      = p.idleTime_s + p.rampEndTime_s;

    % --- Sampling / ADC ---
    p.numADCSamples  = numAdc;
    p.Fs_Hz          = digKsps * 1e3;

    % --- Bandwidths & range resolution ---
    p.B_total_Hz     = p.Slope_Hz_per_s * nanzero(p.rampEndTime_s);      % GUI BW
    if all(~isnan([p.numADCSamples, p.Fs_Hz])) && p.Fs_Hz>0
        T_adc = p.numADCSamples / p.Fs_Hz;
        p.B_sampled_Hz = p.Slope_Hz_per_s * T_adc;                        % FFT BW
        p.rangeRes_m   = 3e8 / (2 * p.B_sampled_Hz);
    else
        p.B_sampled_Hz = NaN; p.rangeRes_m = NaN;
    end

    % --- Frame / Loops / Periodicity ---
    if ~isempty(Raw.AdvancedFrameConfig)
        [numFrames, loops, framePeriod_s] = parse_adv_frame(Raw.AdvancedFrameConfig);
    else
        [numFrames, loops, framePeriod_s] = parse_frame(Raw.FrameConfig);
    end
    % Fallback: parse text labels anywhere in the log if needed
    if isnan(loops) || isnan(numFrames) || isnan(framePeriod_s)
        [loops2, frames2, per2] = parse_frame_from_labels(whole);
        if isnan(loops),        loops = loops2;        end
        if isnan(numFrames),    numFrames = frames2;   end
        if isnan(framePeriod_s),framePeriod_s = per2;  end
    end
    % last-ditch defaults
    if isnan(loops), loops = 1; end
    if isnan(numFrames), numFrames = 1; end

    p.numFrames      = numFrames;
    p.chirpsPerLoop  = loops;
    p.framePeriod_s  = framePeriod_s;

    % --- Chirp segments & chirps per frame (robust) ---
    [segments, nChirpsOnce] = summarize_chirps_robust(Raw.ChirpConfig_ALL, pc);
    p.chirpSegments  = segments;                         % [start end txMask]
    p.chirpsPerFrame = max(1, nChirpsOnce) * max(1, round(loops));

    % --- extras ---
    p.fc_const_raw        = startFreqConst;
    p.hsiClk_Hz_fromLine  = hsiHz;
    p.tsw_sampRate_val    = tswSamp;
    p.lambda_m            = 3e8/77e9;
    p.units_profile_times = u_prof;
    p.notes = 'Args-only parse; Frame labels fallback; robust TX mask & chirp count.';

    RadarParams = p;
end

% ============ helpers ============
function L = read_lines(f)
    raw = fileread(f);
    raw = strrep(raw, sprintf('\r\n'), sprintf('\n'));
    raw = strrep(raw, sprintf('\r'),   sprintf('\n'));
    L = regexp(raw, '\n', 'split')';
    L = L(~cellfun(@isempty,L));
end

function v = get_args_one(L, key)
    v = [];
    idx = find(contains(L, key));
    if isempty(idx), return; end
    S = L{idx(end)};
    m = regexp(S,[regexptranslate('escape',key) '\s*\(([^)]*)\)'],'tokens','once');
    if ~isempty(m), v = parse_args(m{1}); return; end
    m = regexp(S,['API:\s*' regexptranslate('escape',key) '\s*,(.*)$'],'tokens','once');
    if ~isempty(m), v = parse_args(m{1}); return; end
end

function A = get_args_all(L, key)
    A = {};
    idx = find(contains(L, key));
    for ii = 1:numel(idx)
        S = L{idx(ii)};
        m = regexp(S,[regexptranslate('escape',key) '\s*\(([^)]*)\)'],'tokens','once');
        if isempty(m)
            m = regexp(S,['API:\s*' regexptranslate('escape',key) '\s*,(.*)$'],'tokens','once');
        end
        if ~isempty(m), A{end+1} = parse_args(m{1}); end %#ok<AGROW>
    end
end

function nums = parse_args(argstr)
    toks = regexp(argstr, ',', 'split'); toks = strtrim(toks);
    toks = toks(~cellfun(@isempty,toks));
    nums = zeros(1,numel(toks));
    for i=1:numel(toks)
        t = regexprep(toks{i}, '[;)]$', '');
        if ~isempty(regexp(t, '^0x[0-9A-Fa-f]+$', 'once'))
            nums(i) = double(hex2dec(t(3:end)));
        elseif ~isempty(regexp(t, '^[\+\-]?\d+(\.\d+)?([eE][\+\-]?\d+)?$', 'once'))
            nums(i) = str2double(t);
        else
            nums(i) = NaN;
        end
    end
    while ~isempty(nums) && isnan(nums(end)), nums(end)=[]; end
end

function out = first_nonempty(list)
    out = [];
    for i=1:numel(list)
        if ~isempty(list{i}), out = list{i}; return; end
    end
end

function [hsiHz, tswSamp] = parse_tsw(L)
    hsiHz = NaN; tswSamp = NaN;
    idx = find(contains(L, 'Sampling rate'));
    if isempty(idx), return; end
    S = L{idx(end)}; nums = regexp(S, '\d+', 'match');
    if numel(nums)>=2
        hsiHz  = str2double(nums{1});
        tswSamp= str2double(nums{2});
    end
end

function [txMask, rxMask] = interpret_masks(v)
    txMask = uint32(0); rxMask = uint32(0);
    if isempty(v) || numel(v)<2, return; end
    a = sanitize(v(1)); b = sanitize(v(2));
    if bitcount_u32(a)<=3 && bitcount_u32(b)<=4
        txMask=a; rxMask=b;
    elseif bitcount_u32(b)<=3 && bitcount_u32(a)<=4
        txMask=b; rxMask=a;
    else
        if a<=b, txMask=a; rxMask=b; else, txMask=b; rxMask=a; end
    end
end
function u = sanitize(x)
    if isnan(x), u=uint32(0); return; end
    x = floor(abs(x)); u = uint32(mod(x,2^32));
end
function n = bitcount_u32(x)
    n = sum(dec2bin(double(x))=='1');
end
function y = geti(v,i), if numel(v)>=i, y=v(i); else, y=NaN; end, end
function unit = detect_profile_unit(vals)
    vals = vals(~isnan(vals));
    if isempty(vals), unit='us'; return; end
    unit = tern(median(vals) >= 2000, 'tick10ns', 'us');
end
function t = conv_time(x, unit)
    if isnan(x), t=NaN; return; end
    switch unit
        case 'us',       t = x * 1e-6;
        case 'tick10ns', t = x * 10e-9;
        otherwise,       t = NaN;
    end
end
function z = nanzero(x), if isnan(x) || x==0, z=eps; else, z=x; end, end
function s = tern(c,a,b), if c, s=a; else, s=b; end, end

% ---- Frame parsers ----
function [numFrames, loops, framePeriod_s] = parse_frame(fc)
% ar1.FrameConfig(startChirp, endChirp, numLoops, numFrames, periodicity_ms, ...)
    numFrames = NaN; loops = NaN; framePeriod_s = NaN;
    if isempty(fc) || numel(fc)<5, return; end
    loops         = fc(3);
    numFrames     = fc(4);
    framePeriod_s = fc(5) * 1e-3;
end

function [numFrames, loops, framePeriod_s] = parse_adv_frame(af)
    numFrames = NaN; loops = NaN; framePeriod_s = NaN;
    if isempty(af), return; end
    ms  = af(af>=1 & af<=2000);
    tks = af(af>=5e5 & af<=2e8);
    if ~isempty(ms),  framePeriod_s = ms(1)*1e-3; end
    if isempty(ms) && ~isempty(tks), framePeriod_s = tks(1)*10e-9; end
    ints = af(af>0 & af<1e7 & abs(af-round(af))<1e-6);
    if ~isempty(ints), numFrames = ints(1); end
end

function [loops, frames, per_s] = parse_frame_from_labels(text)
% Fallback: read "No of Chirp Loops", "No of Frames", "Periodicity (ms)" anywhere in the log
    loops = NaN; frames = NaN; per_s = NaN;
    m = regexp(text, 'No of Chirp Loops[^0-9]*(\d+)', 'tokens', 'once'); if ~isempty(m), loops  = str2double(m{1}); end
    m = regexp(text, 'No of Frames[^0-9]*(\d+)',      'tokens', 'once'); if ~isempty(m), frames = str2double(m{1}); end
    m = regexp(text, 'Periodicity\s*\(ms\)[^0-9]*(\d+(\.\d+)?)', 'tokens', 'once');
    if ~isempty(m), per_s = str2double(m{1}) * 1e-3; end
end

% ---- Chirp summarizer (robust TX enables) ----
function [segments, nChirpsOnce] = summarize_chirps_robust(ChList, pc)
% Each ar1.ChirpConfig is typically:
% [profId, startCh, endCh, startFreqVar, slopeVar, idleVar, adcStartVar, TX0, TX1, TX2]
% But tail may vary; we search for the *last three* boolean flags after index 7.
    segments = [];
    nChirpsOnce = 0;
    if isempty(ChList), return; end
    profId = pc(1);
    for k=1:numel(ChList)
        v = ChList{k};
        if numel(v) < 3, continue; end
        if v(1) ~= profId, continue; end
        startCh = max(0, round(v(2)));
        endCh   = max(startCh, round(v(3)));              % guard: ensure ≥ start
        txMask  = infer_tx_mask_from_tail(v);
        segments = [segments; startCh, endCh, txMask];    %#ok<AGROW>
        nChirpsOnce = nChirpsOnce + (endCh - startCh + 1);
    end
end

function txMask = infer_tx_mask_from_tail(v)
% find last 3 boolean-ish (0/1) entries after index 7, else fallback to 0/1
    txMask = 0;
    i0 = min(8, numel(v)); tail = v(i0:end);
    isBool = @(x) (~isnan(x)) && (x==0 || x==1);
    bins = [];
    for i=1:numel(tail)
        if isBool(tail(i)), bins = [bins tail(i)]; %#ok<AGROW>
    end
    % keep only the last 3 flags
    if numel(bins) >= 3
        bins = bins(end-2:end);
    elseif numel(tail) >= 3
        % as a fallback, examine the last 3 entries of the vector
        bins = arrayfun(@(x) double(x~=0), tail(max(1,end-2):end));
    else
        bins = [0 0 0];
    end
    % bins assumed ordered [TX0 TX1 TX2]
    for b=1:min(3,numel(bins))
        if bins(b), txMask = bitor(txMask, bitshift(1, b-1)); end
    end
end
