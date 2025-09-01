function cal = calibrate_with_corner(adcData, p, R0_m)
% adcData: [Rx, Ns, Nd, Nf] or [Rx, Ns, Nd] averaged
% p: struct with Slope_Hz_per_s, Fs_Hz, lambda_m, numRx, nfft sizes, etc.
% R0_m: known range to corner reflector

% 1) Average frames (if present)
if ndims(adcData) == 4
    adcData = mean(adcData, 4);
end

% 2) Range-Doppler per-Rx (your own implementation or call helper)
[RD, rAxis, dAxis] = rangeDoppler_perRx(adcData, p); % RD: [Rx, Nr, Nd]

% 3) Find reflector peak near expected range & Doppler~0
[~, r0] = min(abs(rAxis - R0_m));
[~, d0] = min(abs(dAxis - 0));
win = max(1, r0-1):min(numel(rAxis), r0+1);

mag_sum = squeeze(sum(abs(RD(:, win, d0)).^2, 1));
[~, k]   = max(mag_sum);
rHat_m   = rAxis(win(k));

% 4) Range bias
cal.rangeBias_m = rHat_m - R0_m;

% 5) RX relative gain/phase (choose RX1 as ref)
X = squeeze(RD(:, win(k), d0));  % [Rx]
xref = X(1);
cal.rxPhase_rad = angle(X ./ xref);
cal.rxGain_lin  = abs(xref) ./ abs(X);

% 6) TX phase (TDM-MIMO): user should pass adcData split by TX
% Here only placeholder; see docs/primer for extraction by TX burst.
cal.txPhase_rad = [];   % fill if you separate TX bursts

% 7) Optional: boresight angle bias (if AoA module available)
cal.angBias_deg = 0;

% 8) Metadata for reproducibility
cal.meta = p;
cal.meta.R0_m = R0_m;
cal.meta.timestamp = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
cal.meta.hash = DataHash(cal); % use File Exchange DataHash if installed
end
