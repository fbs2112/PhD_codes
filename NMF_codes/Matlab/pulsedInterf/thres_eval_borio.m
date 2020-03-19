clear;
clc;
close all;

addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
load sim_params_2.mat;

alpha = 4.5e11;
deltaT = 12e-6;
params.fs = paramsSignal.Freqsamp;

t = -10e-6 - 1/params.fs:1/params.fs:20e-6 - 1/params.fs;
gaussFun = @(alpha, deltaT, t) exp((-alpha/2).*t.^2) + exp((-alpha/2).*(t - deltaT).^2);
dme = 3*gaussFun(alpha, deltaT, t).';

numberOfRawSamples = floor(paramsSignal.Freqsamp*paramsSignal.Intetime);
delay = 100e-6;
totalSamples = numberOfRawSamples;
numberOfZeros = numberOfRawSamples - length(dme);
monteCarloLoops = 100;
SNR = -25;

params.JNRVector = -5;
paramsSignal.numberOfGPSSignals = 1;

nbits = 4;

for loopIndex = 1:monteCarloLoops
    paramsSignal.FreqDopp = 0;
    interferenceSignal = [zeros(numberOfZeros/4, 1);dme;zeros(numberOfZeros*3/4, 1)];
    Timeofthisloop = 0:totalSamples-1;
    Carrphase = mod(2*pi*(paramsSignal.FreqDopp)*Timeofthisloop/paramsSignal.Freqsamp,2*pi);
    Carrier = exp(1i*Carrphase).';
    
    interferenceSignal = interferenceSignal.*Carrier;
    paramsSignal.FreqDopp = 0;
    GPSSignals = GPSGen(paramsSignal);
    GPSSignals = GPSSignals(1:numberOfRawSamples,:);
    
    GPSSignals = [zeros(round(delay*params.fs), 1);GPSSignals];
    GPSSignalsPower = pow_eval(GPSSignals);
    interferenceSignal = [zeros(round(delay*params.fs), 1);interferenceSignal];
    interferenceSignalPower = pow_eval(interferenceSignal);
%     noise(:,loopIndex) = randn(length(GPSSignals), 1) + 1j*randn(length(GPSSignals), 1);
    noise(:,loopIndex) = randn(length(GPSSignals), 1);
    noisePower = pow_eval(noise(:,loopIndex));
    varN(loopIndex) = var(noise(:,loopIndex));

    GPSSignalsAux = GPSSignals;
    interferenceSignalAux = interferenceSignal;
    GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
    mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise(:,loopIndex);
    mixtureGPSQuant(:,loopIndex) = quantise_gps(mixtureGPS, nbits, 1);
    interferenceSignalAux = interferenceSignalAux*sqrt(noisePower*10^(params.JNRVector/10)/interferenceSignalPower);
    mixtureSignal(:,loopIndex) = mixtureGPS + interferenceSignalAux;
    mixtureSignalQuant(:,loopIndex) = quantise_gps(mixtureSignal(:,loopIndex), nbits, 1);
end

histogram(abs(noise(:)));
figure;
histogram(abs(mixtureSignal(:)));
figure;
histogram(abs(mixtureSignalQuant(:)));

figure;
histogram(abs(mixtureGPSQuant(:)));

% varN = 2;
% Pfa = 0.01;
% threshold = sqrt(-2*varN*log(Pfa));

rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);