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
    nav(i,1) = (PRN.*GPSMultiplier)*(PRN.*GPSMultiplier).';

    totalSamples = length(mixtureGPS);
    Timeofthisloop = 0:totalSamples-1;
    Carrphase = mod(2*pi*(paramsSignal.FreqDopp)*Timeofthisloop/paramsSignal.Freqsamp,2*pi);
    Carrier = exp(1i*Carrphase).';

    mixtureGPSNoDoppler = (mixtureGPS.*conj(Carrier));

    mixtureGPSNav(i,1) = mixtureGPSNoDoppler.'*PRN.';
   
end
trueSNR = pow2db(pow_eval(nav) / (noisePower*monteCarloLoops));

signalPower = mixtureGPSNav'*mixtureGPSNav/length(mixtureGPSNav);
gpsSignalPowerHat = mean(abs(real(mixtureGPSNav)))^2;

SNR_SNV = pow2db((gpsSignalPowerHat/(signalPower - gpsSignalPowerHat)));

Ta = Tc;
M = length(mixtureGPS);

Pn = sum(real(mixtureGPSNav))^2 + sum(imag(mixtureGPSNav))^2;
Pw = sum(real(mixtureGPSNav).^2 + imag(mixtureGPSNav).^2);

SNR_NWPR = pow2db((M*(Pn/Pw) - 1) / ((M - (Pn/Pw))));

M4 = mean(abs(mixtureGPSNav).^4);
Pd = sqrt(2*(signalPower^2) - M4);
SNR_MM = pow2db((Pd/(signalPower - Pd)));

rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);