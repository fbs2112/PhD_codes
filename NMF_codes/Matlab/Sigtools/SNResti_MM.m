function SNRestimate = SNResti_MM(Signal)

signalPower = Signal'*Signal/length(Signal);

M4 = mean(abs(Signal).^4);
Pd = sqrt(2*(signalPower^2) - M4);
SNRestimate = pow2db((Pd/(signalPower - Pd)));


