clear;
clc;
close all;

addpath(['..' filesep 'Sigtools' filesep])
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);

load sim_params_2.mat;

SNR = -25;
Tc = paramsSignal.Intetime;
monteCarloLoops = 100;
paramsSignal.numberOfGPSSignals = 1;
[GPSSignals, PRN] = GPSGen(paramsSignal);
GPSSignalsPower = pow_eval(GPSSignals);

for i = 1:monteCarloLoops
    noise = randn(length(GPSSignals), 1) + 1j*randn(length(GPSSignals), 1);
    noisePower = pow_eval(noise);

    GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
    mixtureGPS = sum(GPSSignals.*GPSMultiplier, 2) + noise;

    trueCN0(i) = pow2db((pow_eval((PRN.*PRN).')/noisePower)/Tc);

    totalSamples = length(mixtureGPS);
    Timeofthisloop = 0:totalSamples-1;
    Carrphase = mod(2*pi*(paramsSignal.FreqDopp)*Timeofthisloop/paramsSignal.Freqsamp,2*pi);
    Carrier = exp(1i*Carrphase).';

    mixtureGPSNoDoppler = (mixtureGPS.*conj(Carrier));

    mixtureGPSNav = mixtureGPSNoDoppler.*PRN.';
    signalPower = mixtureGPSNav'*mixtureGPSNav/length(mixtureGPSNav);
    gpsSignalPowerHat = mean(abs(real(mixtureGPSNav)))^2;

    CN0_SNV(i) = pow2db((gpsSignalPowerHat/(signalPower - gpsSignalPowerHat))/Tc);

    Ta = Tc;
    M = length(mixtureGPS);

    Pn = sum(real(mixtureGPSNav))^2 + sum(imag(mixtureGPSNav))^2;
    Pw = sum(real(mixtureGPSNav).^2 + imag(mixtureGPSNav).^2);

    CN0_NWPR(i) = pow2db((M*(Pn/Pw) - 1) / (Ta * (M - (Pn/Pw))));
    % M4 = mean(abs(mixtureGPSNav).^4);
    % Pd = sqrt(2*(signalPower^2) - M4);
    % CN0_MM = pow2db((Pd/(signalPower - Pd))/Tc);
end
CN0_SNV_av = mean(CN0_SNV);
CN0_NWPR_av = mean(CN0_NWPR);
trueCN0_av = mean(trueCN0); 

rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);