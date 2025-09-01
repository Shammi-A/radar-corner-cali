function Xc = apply_rx_tx_cal(X, cal)
% X: [Rx, ...] complex data in any later stage (range/Doppler/AoA)
G = reshape(cal.rxGain_lin, [], 1);
P = reshape(exp(-1j*cal.rxPhase_rad), [], 1);
Xc = bsxfun(@times, bsxfun(@times, X, G), P);
end
