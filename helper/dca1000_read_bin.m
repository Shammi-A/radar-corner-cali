function [adc, meta] = dca1000_read_bin(binPath, p)
% dca1000_read_bin  (Drop-in shim) Read AWR2243+DCA1000 BIN using Studio "Format 0".
% Usage:
%   [adc, meta] = dca1000_read_bin(binPath, p)
% where p is your parsed config struct (from your existing demo_parse_log).

% ---- Map fields from p -> cfg (robust + nested + inference) ----
cfg.numAdcSamples = pickDeep(p, {'numADCSamples','ADC_samples','numSamples','adcSamples','numAdcSamples'});
cfg.numRx         = pickDeep(p, {'numRx','numRX','rx','numReceiveAntennas','numReceivers'});
cfg.numTx         = pickDeep(p, {'numTx','tx','numTransmitAntennas','numTxAntennas','numTransmitters'});

% Loops: try explicit first; else infer (with I/Q double-count guard)
cfg.numLoops = pickDeep(p, {'numLoops','loops','numChirpLoops','numChirpsLoop','numLoopsPerFrame','chirpLoops','loopCount'});
if isempty(cfg.numLoops)
    cpf = pickDeep(p, {'numChirpsPerFrame','chirpsPerFrame','M','chirps_per_frame','ChirpsPerFrame'});
    if ~isempty(cpf) && ~isempty(cfg.numTx)
        cpf = double(cpf); ntx = double(cfg.numTx);
        if mod(cpf, 2*ntx) == 0 && cpf/(2*ntx) <= 4096   % guard typical ranges
            L = cpf/(2*ntx);   % correct common 768 (I/Q double count) -> 128
        else
            L = cpf/ntx;
        end
        if abs(L - round(L)) < 1e-6, cfg.numLoops = round(L); end
        fprintf('ℹ️  Inferred numLoops = %d (from chirpsPerFrame=%d, numTx=%d)\n', cfg.numLoops, cpf, cfg.numTx);
    end
end

% Final assertions
assert(~isempty(cfg.numAdcSamples), 'Missing ADC samples in p (e.g., p.numADCSamples)');
assert(~isempty(cfg.numRx),        'Missing numRx in p (e.g., p.numRx)');
assert(~isempty(cfg.numTx),        'Missing numTx in p (e.g., p.numTx)');
assert(~isempty(cfg.numLoops),     'Missing numLoops/loops in p. Add p.numLoops=128 or p.numChirpsPerFrame.');

cfg.verbose = true;

% ---- Core read (Format 0: 16-bit I + 16-bit Q) ----
[adc, meta] = fmt0_reader(binPath, cfg);

% ---- Optional: frame sanity print ----
expectedCpf = cfg.numLoops * cfg.numTx;  % e.g., 128*3=384
if abs(meta.numFrames_est - round(meta.numFrames_est)) < 1e-3
    fprintf('Frames (est): %d\n', round(meta.numFrames_est));
else
    warning('Not an integer number of frames (est = %.3f). Check Chirp/Frame setup.', meta.numFrames_est);
end
end

% =======================================================================
function v = pickDeep(s, names)
% pickDeep: search top-level, then 1-level nested structs for first hit.
v = [];
if ~isstruct(s), return; end
for i = 1:numel(names)
    n = names{i};
    if isfield(s, n) && ~isempty(s.(n)), v = s.(n); return; end
    fns = fieldnames(s);
    for k = 1:numel(fns)
        sub = s.(fns{k});
        if isstruct(sub) && isfield(sub, n) && ~isempty(sub.(n))
            v = sub.(n); return;
        end
    end
end
end

% =======================================================================
function [adc, meta] = fmt0_reader(binPath, cfg)
% fmt0_reader  Read BIN in mmWave Studio "Format 0" (16-bit I/Q).

if ~isfield(cfg,'verbose'), cfg.verbose = true; end
Nr  = double(cfg.numAdcSamples);
Nrx = double(cfg.numRx);
Ntx = double(cfg.numTx);
Nlp = double(cfg.numLoops);
vb  = cfg.verbose;

bytesPerInt16 = 2;
comps         = 2;  % I+Q
bytesPerSampPerRx = comps * bytesPerInt16;              % 4 bytes
bytesPerChirp     = Nr * Nrx * bytesPerSampPerRx;       % e.g., 256*4*4 = 4096
chirpsPerFrame_expected = Nlp * Ntx;                    % e.g., 384
bytesPerFrame     = chirpsPerFrame_expected * bytesPerChirp;

% -- Read raw bytes
fid = fopen(binPath, 'rb'); assert(fid>0, 'Could not open %s', binPath);
raw = fread(fid, inf, 'uint8=>uint8'); fclose(fid);
fileBytes = numel(raw);

if vb
    fprintf('📥 %s | %d bytes\n', binPath, fileBytes);
    fprintf('   Bytes/chirp=%d | Expected chirps/frame=%d | Bytes/frame=%d\n', ...
        bytesPerChirp, chirpsPerFrame_expected, bytesPerFrame);
end

% -- Trim to whole chirps (drop dangling tail safely)
remBytes_chirp = mod(fileBytes, bytesPerChirp);
if remBytes_chirp ~= 0
    if vb, fprintf('⚠️  Tail not a multiple of chirp (%d bytes). Dropping tail.\n', remBytes_chirp); end
    raw = raw(1:end-remBytes_chirp);
    fileBytes = numel(raw);
end

% -- Convert to int16 stream and reshape
int16stream = typecast(raw, 'int16'); clear raw;
int16sPerChirp = Nr * (2*Nrx);
totalChirps    = numel(int16stream) / int16sPerChirp;
assert(abs(totalChirps - floor(totalChirps)) < 1e-9, 'Internal size error.');
totalChirps    = floor(totalChirps);

x = reshape(int16stream, [2*Nrx, Nr, totalChirps]); clear int16stream;
I = x(1:2:end, :, :); Q = x(2:2:end, :, :);
adc = complex(single(I), single(Q)); clear I Q x;

% -- Frame-alignment sanity
remBytes_frame = mod(fileBytes, bytesPerFrame);
if vb
    if remBytes_frame == 0
        fprintf('✅ File aligns to whole frames (%d bytes/frame).\n', bytesPerFrame);
    else
        fprintf('ℹ️  File not a whole number of frames; remainder %d bytes (OK to proceed).\n', remBytes_frame);
    end
end

% -- Meta
meta = struct();
meta.fileBytes               = fileBytes;
meta.bytesPerChirp           = bytesPerChirp;
meta.chirpsPerFrame_expected = chirpsPerFrame_expected;
meta.bytesPerFrame           = bytesPerFrame;
meta.numChirps               = totalChirps;
meta.numFrames_est           = fileBytes / bytesPerFrame; % may be fractional
meta.tailDropped_bytes       = remBytes_chirp;
meta.format                  = 'Format0_I16Q16';
meta.shape_adc               = size(adc); % [Nrx Nr Nchirps]
end
