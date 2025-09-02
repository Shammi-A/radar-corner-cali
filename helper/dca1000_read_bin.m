function out = dca1000_read_bin(binFile, RadarParams, Raw, opts)
% DCA1000 reader using parsed RadarParams/Raw from mmws_parse_log.m
% Output:
%   out.cube      : [numRx x numADCSamples x chirpsPerFrame x numFramesUsed] (complex or real)
%   out.meta      : struct with parsed/derived fields
%
% Usage:
%   out = dca1000_read_bin('adc.bin', outParse.RadarParams, outParse.Raw);

    arguments
        binFile (1,:) char
        RadarParams struct
        Raw struct = struct()
        opts.forceComplex (1,1) logical = false   % force IQ mode
        opts.forceReal    (1,1) logical = false   % force real mode
        opts.returnClass  (1,:) char   = 'double' % 'double' | 'single' | 'int16'
    end

    assert(exist(binFile,'file')==2, 'Bin file not found: %s', binFile);

    % ---- pull key config from parsed params ----
    numRx    = RadarParams.numRx;
    Ns       = RadarParams.numADCSamples;
    Nc_pf    = RadarParams.chirpsPerFrame;        % chirps per frame (across TX pattern × loops)
    Nf_conf  = RadarParams.numFrames;             % configured frames (may be > actually captured)

    assert(numRx>0 && Ns>0 && Nc_pf>0, 'Invalid parsed config (numRx=%d, Ns=%d, chirpsPerFrame=%d).', numRx, Ns, Nc_pf);

    % ---- derive ADC bits & guess complex/real from Raw.* (with a fallback) ----
    adcBits = 16;    % default; override if Raw.AdcOutConfig says otherwise
    isComplex_hint = []; % [] = unknown; true/false if we can infer

    if isfield(Raw,'AdcOutConfig') && ~isempty(Raw.AdcOutConfig)
        v = Raw.AdcOutConfig;
        % Common conventions in Studio logs:
        %  v(2): ADC bits (0=12, 1=14, 2=16)  -> your logs show '2' => 16-bit
        %  v(1): format (0=real, 1/2=complex) varies across DFPs; treat 2 as complex if present.
        if numel(v) >= 2 && ~isnan(v(2)), adcBits = pick_bits(v(2)); end
        if numel(v) >= 1 && ~isnan(v(1))
            if v(1) >= 1, isComplex_hint = true; else, isComplex_hint = false; end
        end
    end

    % ---- read raw int16 (little-endian as produced by DCA1000) ----
    fid = fopen(binFile, 'r');
    raw = fread(fid, inf, 'int16=>int16', 0, 'ieee-le');
    fclose(fid);

    if adcBits ~= 16
        raw = sign_extend_int16(raw, adcBits);   % handle 12/14-bit packing
    end

    Nint16 = numel(raw);

    % ---- choose complex vs real (auto unless forced) ----
    if opts.forceComplex && opts.forceReal
        error('Choose either forceComplex or forceReal, not both.');
    end
    if opts.forceComplex
        isComplex = true;
    elseif opts.forceReal
        isComplex = false;
    else
        isComplex = autodetect_cplx(Nint16, numRx, Ns, Nc_pf, isComplex_hint);
    end

    wordsPerSample = 2 * double(isComplex) + 1 * double(~isComplex);  % 2 int16 per sample if complex

    % ---- infer #frames actually captured from file size ----
    wordsPerFrame = numRx * Ns * Nc_pf * wordsPerSample;
    if wordsPerFrame == 0
        error('wordsPerFrame computed as zero (check config).');
    end
    Nf_file = floor(double(Nint16) / double(wordsPerFrame));
    remWords = double(Nint16) - Nf_file*wordsPerFrame;

    if remWords ~= 0
        warning('Bin size not a multiple of one frame: dropping %d trailing int16.', remWords);
    end

    % ---- practical #frames to load ----
    if isfinite(Nf_conf) && ~isnan(Nf_conf) && Nf_conf > 0
        Nf = min(Nf_conf, Nf_file);
    else
        Nf = Nf_file;
    end
    if Nf == 0
        error('No full frames present in file with current config (Nint16=%d, perFrame=%d).', Nint16, wordsPerFrame);
    end

    % ---- trim and vectorize ----
    raw = raw(1 : Nf * wordsPerFrame);

    if isComplex
        % TI Studio packing (common): I0, I1, Q0, Q1, I2, I3, Q2, Q3, ...
        % Turn raw int16 -> complex vector of length (Nint16/2)
        v = complex( zeros(numel(raw)/2,1) );
        v(1:2:end) = complex(double(raw(1:4:end)), double(raw(3:4:end)));
        v(2:2:end) = complex(double(raw(2:4:end)), double(raw(4:4:end)));

        % reshape to [Ns*numRx, totalChirps] then to [numRx, Ns, Nc_pf, Nf]
        totalChirps = Nc_pf * Nf;
        if numel(v) ~= Ns*numRx*totalChirps
            error('Complex length mismatch. Got %d complex samples, expect %d.', numel(v), Ns*numRx*totalChirps);
        end
        M = reshape(v, Ns*numRx, totalChirps);
        C = reshape(permute(reshape(M, Ns, numRx, totalChirps), [2 1 3]), numRx, Ns, Nc_pf, Nf);

        cube = cast(C, opts.returnClass);

    else
        % Real-only: one int16 per ADC sample
        v = double(raw);

        totalChirps = Nc_pf * Nf;
        if numel(v) ~= Ns*numRx*totalChirps
            error('Real length mismatch. Got %d samples, expect %d.', numel(v), Ns*numRx*totalChirps);
        end

        M = reshape(v, Ns*numRx, totalChirps);
        C = reshape(permute(reshape(M, Ns, numRx, totalChirps), [2 1 3]), numRx, Ns, Nc_pf, Nf);

        % keep real as requested class
        cube = cast(C, opts.returnClass);
    end

    % ---- assemble output ----
    out = struct();
    out.cube = cube;
    out.meta = struct( ...
        'binFile',        binFile, ...
        'isComplex',      isComplex, ...
        'adcBits',        adcBits, ...
        'numRx',          numRx, ...
        'numADCSamples',  Ns, ...
        'chirpsPerFrame', Nc_pf, ...
        'framesInFile',   Nf_file, ...
        'framesUsed',     Nf, ...
        'Fs_Hz',          RadarParams.Fs_Hz, ...
        'slope_MHz_per_us', RadarParams.slope_MHz_per_us, ...
        'framePeriod_s',  RadarParams.framePeriod_s, ...
        'frameRate_Hz',   RadarParams.frameRate_Hz ...
    );
end

% ---------- helpers ----------
function bits = pick_bits(code)
% Map common Studio codes to ADC bits (best-effort)
    switch round(code)
        case 0, bits = 12;
        case 1, bits = 14;
        otherwise, bits = 16;  % 2 or anything else → 16-bit
    end
end

function raw16 = sign_extend_int16(raw16, validBits)
% Sign-extend raw16 that contain only "validBits" LSBs of signal
% (typical when ADC is 12/14-bit captured into 16-bit containers)
    mask    = int16(bitshift(int16(1), validBits) - 1);
    signbit = int16(bitshift(int16(1), validBits-1));
    x = bitand(raw16, mask);
    neg = bitand(x, signbit) ~= 0;
    x(neg) = x(neg) - bitshift(int16(1), validBits);
    raw16 = x;
end

function isCx = autodetect_cplx(Nint16, numRx, Ns, Nc_pf, hint)
% Decide complex vs real based on file length & optional hint
    if ~isempty(hint)
        isCx = logical(hint);
        return;
    end
    wordsPerFrame_real = numRx * Ns * Nc_pf;
    wordsPerFrame_cx   = wordsPerFrame_real * 2;

    r_real = mod(Nint16, wordsPerFrame_real);
    r_cx   = mod(Nint16, wordsPerFrame_cx);

    if r_cx == 0 && r_real ~= 0
        isCx = true;
    elseif r_real == 0 && r_cx ~= 0
        isCx = false;
    elseif r_cx == 0 && r_real == 0
        % both divide cleanly — prefer complex (most common with LVDS IQ)
        isCx = true;
    else
        % neither divides cleanly: choose the smaller remainder (best fit)
        isCx = (r_cx <= r_real);
    end
end
