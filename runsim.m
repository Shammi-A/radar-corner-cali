clc; close all; clear;
rehash;
%% select the test
test = "T4"; % options:( T1, T2, T3, T4)

outParse = demo_parse_log("G:\My Drive\mmWave -sensing-data\calibration", test);
p   = outParse.RadarParams;
raw = outParse.Raw;
binPath = outParse.binFiles.path{1};     % assumes your helper returns a table with 'path'

%% Read the DCA1000 .bin using the parsed config
r = dca1000_read_bin(binPath, p, raw);

% r.cube is [numRx x Ns x chirpsPerFrame x frames] complex (or real)
rx1_frame1 = squeeze(r.cube(1,:,:,1));   % [Ns x chirpsPerFrame]

%% Calibration --


%% Corner distance (meters) —ground-truth here

% Value in feet
feet = 12; % options( 3, 6, 9, 12)

% Conversion factor
ft2m = 0.3048;

Rcorner =  feet * ft2m; % change as i move the ground-truth 

% Calibrate
cal = corner_calibrate(r.cube, p, 'cornerRange_m', Rcorner, ...
                                'rangeSearch_m', Rcorner + [-0.6 0.6], ...
                                'plot', true);

% Apply calibration to the same cube (or save 'cal' and apply to future captures)
cubeCal = apply_corner_calib(r.cube, p, cal);

% After you have r.cube, p, and cal from corner_calibrate:
stats = corner_calib_check(r.cube, p, cal, ...
    'useFrames', 40, ...
    'cornerRange_m', Rcorner, ...   
    'plot', true);

% Inspect numeric results
disp(stats);
