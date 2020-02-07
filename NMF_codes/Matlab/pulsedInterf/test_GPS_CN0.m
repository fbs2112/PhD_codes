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
    nav(i,1) = (PRN.*GPSMultiplier)*(PRN).';

    totalSamples = length(mixtureGPS);
    Timeofthisloop = 0:totalSamples-1;
    Carrphase = mod(2*pi*(paramsSignal.FreqDopp)*Timeofthisloop/paramsSignal.Freqsamp,2*pi);
    Carrier = exp(1i*Carrphase).';

    mixtureGPSNoDoppler = (mixtureGPS.*conj(Carrier));

    mixtureGPSNav(i,1) = mixtureGPSNoDoppler.'*PRN.';
   
end
trueSNR = pow2db(pow_eval(nav) / (noisePower*monteCarloLoops.^2));

signalPower = mixtureGPSNav'*mixtureGPSNav/length(mixtureGPSNav);
gpsSignalPowerHat = mean(abs(real(mixtureGPSNav)))^2;

SNR_SNV = pow2db((gpsSignalPowerHat/(signalPower - gpsSignalPowerHat)));

Ta = Tc;
M = 20;

mixtureGPSNavBuff = buffer(mixtureGPSNav, M);

for i = 1:size(mixtureGPSNavBuff, 2)
    Pn(i) = sum(real(mixtureGPSNavBuff(:,i)))^2 + sum(imag(mixtureGPSNavBuff(:,i)))^2;
    Pw(i) = sum(real(mixtureGPSNavBuff(:,i)).^2 + imag(mixtureGPSNavBuff(:,i)).^2);
end
averagePn = mean(Pn);
averagePw = mean(Pw);

SNR_NWPR = pow2db((M*(averagePn/averagePw) - 1) / ((M - (averagePn/averagePw))));

M4 = mean(abs(mixtureGPSNav).^4);
Pd = sqrt(2*(signalPower^2) - M4);
SNR_MM = pow2db((Pd/(signalPower - Pd)));

SNR_pai = SNResti_pai(mixtureGPSNav);

rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);