function mmws_print_summary(p)
% Human-readable summary of RadarParams (AWR2243 + DCA1000)

    fprintf('--- Parsed Radar Configuration ---\n');

    if isfield(p,'txChanMask')
        fprintf('TX enabled mask: 0x%X  | numTx: %d\n', p.txChanMask, p.numTx);
    end
    if isfield(p,'rxChanMask')
        fprintf('RX enabled mask: 0x%X  | numRx: %d\n', p.rxChanMask, p.numRx);
    end

    if isfield(p,'numADCSamples') && isfield(p,'Fs_Hz')
        fprintf('ADC samples: %d   Fs: %.3f MHz\n', p.numADCSamples, p.Fs_Hz/1e6);
    end

    fprintf('Times (s): idle=%.6g, adcStart=%.6g, rampEnd=%.6g, Tchirp=%.6g  [%s]\n', ...
        p.idleTime_s, p.adcStartTime_s, p.rampEndTime_s, p.T_chirp_s, p.units_profile_times);

    fprintf('Slope: %.3f MHz/us  (%.3e Hz/s)\n', p.slope_MHz_per_us, p.Slope_Hz_per_s);

    fprintf('Sweep (Hz): total=%.3e  sampled=%.3e\n', p.B_total_Hz, p.B_sampled_Hz);
    if isfinite(p.rangeRes_m)
        fprintf('Range resolution: %.3f m (%.1f cm)\n', p.rangeRes_m, 100*p.rangeRes_m);
    else
        fprintf('Range resolution: NaN\n');
    end

    if isfield(p,'numFrames') && isfield(p,'chirpsPerLoop') && isfield(p,'framePeriod_s')
        fprintf('Frame: frames=%d  | loops=%d  | period=%.3f ms  | rate=%.3f Hz\n', ...
            p.numFrames, p.chirpsPerLoop, 1e3*p.framePeriod_s, p.frameRate_Hz);
    end

    if isfield(p,'chirpsPerFrame')
        fprintf('Chirps per frame (computed): %d\n', p.chirpsPerFrame);
    end

    if isfield(p,'dopplerBins') && isfield(p,'vRes_mps') && isfield(p,'vMax_mps')
        fprintf('Doppler: Ndopp=%d  | Δv=%.3f m/s  | v_max=±%.2f m/s\n', ...
            p.dopplerBins, p.vRes_mps, p.vMax_mps);
    end

    if isfield(p,'chirpSegments') && ~isempty(p.chirpSegments)
        seg = p.chirpSegments;
        for i = 1:size(seg,1)
            fprintf('  Chirp segment %d: start=%d end=%d txMask=0x%X\n', i, seg(i,1), seg(i,2), seg(i,3));
        end
    end

    if isfield(p,'hsiClk_Hz_fromLine') && ~isnan(p.hsiClk_Hz_fromLine)
        fprintf('HSI clock (from log text): %.0f MHz\n', p.hsiClk_Hz_fromLine/1e6);
    end
end
