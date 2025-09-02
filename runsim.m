% runsim.m  —  Stream one BIN, reshape to frames, and run corner calibration.
clc; close all; clear;
rehash;

%% ---------------- User inputs ---------------- %%
test   = "T5";          % options: T1, T2, T3, T4, T5
feet   = 15;            % ground-truth corner distance in feet (3, 6, 9, 12, 15, ...)
root   = "G:\My Drive\mmWave -sensing-data\calibration";

%% -------------- Parse log + locate BIN -------------- %%
outParse = demo_parse_log(root, test);   % keep your existing version
p        = outParse.RadarParams;
binPath  = outParse.binFiles.path{1};    % first *_Raw_*.bin

% Inject any fields your parser might have missed (taken from your log printout)
p.numTx          = set_default(p, 'numTx', 3);
p.numRx          = set_default(p, 'numRx', 4);
p.numLoops       = set_default(p, 'numLoops', 128);
p.numADCSamples  = set_default(p, 'numADCSamples', 256);

%% -------------- Read BIN (Format-0: 16-bit I/Q) -------------- %%
[adc, meta] = dca1000_read_bin(binPath, p);
% adc: [numRx x numSamples x numChirpsTotal]

% Sanity: expected chirps/frame from config
cpf = meta.chirpsPerFrame_expected;   % = p.numLoops * p.numTx (e.g., 384)

%% -------------- Align to whole frames & reshape -------------- %%
Nrx = size(adc,1); Nr = size(adc,2); Nc = size(adc,3);

Nc_valid = floor(Nc / cpf) * cpf;     % drop partial frame at tail, if any
if Nc_valid < Nc
    fprintf('Trimming %d chirps to align to whole frames.\n', Nc - Nc_valid);
    adc = adc(:,:,1:Nc_valid);
end

nFrames = Nc_valid / cpf;
cube    = reshape(adc, [Nrx, Nr, cpf, nFrames]);  % [Rx x Ns x chirpsPerFrame x frames]

% Keep old API so your corner_* functions work unchanged
r        = struct();
r.cube   = cube;
r.meta   = meta;
r.p      = p;
% HSI/data link rate (from DataConfig screenshot: 600 Mbps DDR)
if ~isfield(p,'HSI_Hz') || isempty(p.HSI_Hz) || p.HSI_Hz==0
    p.HSI_Hz = 600e6;     % per-lane bit clock; informational only
end

% Optional quick peek (first frame, Rx1)
% rx1_frame1 = squeeze(cube(1,:,:,1));  % [Ns x chirpsPerFrame]

%% -------------- Corner calibration -------------- %%
ft2m    = 0.3048;
Rcorner = feet * ft2m;

% Calibrate (expects r.cube-style input)
cal = corner_calibrate(r.cube, p, ...
        'cornerRange_m', Rcorner, ...
        'rangeSearch_m', Rcorner + [-0.6 0.6], ...
        'plot', true);

% Apply calibration to the same cube (or save 'cal' for future runs)
cubeCal = apply_corner_calib(r.cube, p, cal); %#ok<NASGU>

% Validate calibration quality
stats = corner_calib_check(r.cube, p, cal, ...
            'useFrames', 64, ...
            'cornerRange_m', Rcorner, ...
            'plot', true);

disp(stats);

%% -------------- Helpers -------------- %%
function v = set_default(s, field, fallback)
    if isfield(s, field) && ~isempty(s.(field))
        v = s.(field);
    else
        v = fallback;
        fprintf('ℹ️  Using default %s = %s\n', field, mat2str(v));
    end
end
