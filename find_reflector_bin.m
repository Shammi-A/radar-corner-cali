function [rIdx, dIdx] = find_reflector_bin(RD, rAxis, dAxis, R0)
[~, r0] = min(abs(rAxis - R0));
[~, d0] = min(abs(dAxis - 0));
win = max(1, r0-1):min(numel(rAxis), r0+1);
mag = squeeze(sum(abs(RD(:, win, d0)).^2,1));
[~,k] = max(mag);
rIdx = win(k); dIdx = d0;
end
