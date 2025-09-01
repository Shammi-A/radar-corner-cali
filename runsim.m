clc; close all; clear;
rehash;
addpath('helper');   % <— absolute path so it works from anywhere

rootDir = "G:\My Drive\mmWave -sensing-data\calibration";

out = demo_parse_log(rootDir, 'T1');

% Parsed params:
p = out.RadarParams;
raw = out.Raw;
bins = out.binFiles;
