function SNRestimate = SNResti_SNV(Signal)

signalPower = Signal'*Signal/length(Signal);
gpsSignalPowerHat = mean(abs(real(Signal)))^2;
SNRestimate = pow2db((gpsSignalPowerHat/(signalPower - gpsSignalPowerHat)));

