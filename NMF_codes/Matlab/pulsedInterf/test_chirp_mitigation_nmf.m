clear;
clc;
close all;

addpath(['..' filesep 'Sigtools' filesep])
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
addpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['.' filesep 'data']);

load nmf_training_26.mat;
load sim_params_3.mat;

monteCarloLoops = 1;
SNR = -25;
nbits = 0;

paramsNMF1.JNRVector = [-5 0 10 30 50];
paramsNMF1.JNRVector = 10;

JNRVector = paramsNMF1.JNRVector;
paramsNMF1.fs = paramsSignal.Freqsamp;
paramsNMF1.nfft = 256;
paramsNMF1.nperseg = 256;
paramsNMF1.overlap = paramsNMF1.nperseg - 1;
paramsNMF1.hop_size = paramsNMF1.nperseg - paramsNMF1.overlap;
paramsNMF1.window = ones(paramsNMF1.nperseg, 1);
paramsNMF1.specType = 'power';
paramsNMF1.numberOfSources = 5;
paramsNMF1.init = 'random';
paramsNMF1.betaDivergence = 'kullback-leibler';
paramsNMF1.numberOfIterations = 500;
paramsNMF1.tolChange = 1e-6;
paramsNMF1.tolError = 1e-6;
paramsNMF1.repetitions = 1;
paramsNMF1.verbose = false;
paramsNMF1.transform = true;
paramsNMF1.semi = false;

numberOfRawSamples = floor(paramsSignal.Freqsamp*paramsSignal.Intetime);
delay = 10e-6;
totalSamples = numberOfRawSamples;

bandwidthVector = (2:3:14)*1e6;

periodVector = 8.62e-6;

paramsSignal.Noneperiod = round(periodVector*paramsNMF1.fs);                   % number of samples with a sweep time
paramsSignal.Initphase = 0;

Timeofthisloop = 0:totalSamples-1;
Carrphase = mod(2*pi*(paramsSignal.FreqDopp)*Timeofthisloop/paramsSignal.Freqsamp,2*pi);
Carrier = exp(1i*Carrphase).';

mixtureSignal = zeros(totalSamples, length(paramsNMF1.JNRVector), length(nbits), monteCarloLoops);
xHat = zeros(totalSamples, 2, length(paramsNMF1.JNRVector), length(nbits), length(bandwidthVector), monteCarloLoops);

paramsSignal.FreqDopp = 1e3;
paramsSignal.numberOfGPSSignals = 1;

GPSSignals = GPSGen(paramsSignal);
GPSSignals = GPSSignals(1:numberOfRawSamples,:);
GPSSignals = [GPSSignals(end - round(delay*paramsNMF1.fs)+1:end,:);GPSSignals(1:end - round(delay*paramsNMF1.fs),:)]; % Introducing artificial code delay
GPSSignalsPower = pow_eval(GPSSignals);

paramsNMF2 = paramsNMF1;
paramsNMF2.numberOfSources = paramsNMF1.numberOfSources*2;
paramsNMF2.init = 'custom';
paramsNMF2.transform = false;
paramsNMF2.semi = false;

[PxxGPS, f, t] = spectrogram(GPSSignals, paramsNMF1.window, paramsNMF1.overlap, paramsNMF1.nfft, paramsNMF1.fs, 'centered', paramsNMF1.specType);

for loopIndex = 1:monteCarloLoops
    loopIndex
    
    if paramsSignal.FreqDopp ~= 0
        noise = randn(totalSamples, 1) + 1j*randn(totalSamples, 1);
        noiseVar = 2;
    else
        noise = randn(totalSamples, 1);
        noiseVar = 1;
    end
    noisePower = pow_eval(noise);
    
    for bandwidthIndex = 1:1%length(bandwidthVector)
        bandwidthIndex
        paramsSignal.IFmin = -bandwidthVector(bandwidthIndex)/2;                                                  % start frequency
        paramsSignal.IFmax = bandwidthVector(bandwidthIndex)/2;                    % end frequency
        paramsSignal.foneperiod(1:paramsSignal.Noneperiod) = linspace(paramsSignal.IFmin, paramsSignal.IFmax, paramsSignal.Noneperiod);
        
        [interferenceSignal, ~] = interferenceGen(paramsSignal);
        interferenceSignal = interferenceSignal(1:numberOfRawSamples);
        interferenceSignal = interferenceSignal.*Carrier;
        interferenceSignalPower = pow_eval(interferenceSignal);
        
        for nbitsIndex = 1:length(nbits)
            
            for JNRIndex = 1:length(paramsNMF1.JNRVector)
                GPSSignalsAux = GPSSignals;
                interferenceSignalAux = interferenceSignal;
                GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
                mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
                interferenceSignalAux = interferenceSignalAux*sqrt(noisePower*10^(paramsNMF1.JNRVector(JNRIndex)/10)/interferenceSignalPower);
                mixtureSignal(:,JNRIndex,nbitsIndex,loopIndex) = mixtureGPS + interferenceSignalAux;
            end
           paramsNMF2.W0 = W0{1,bandwidthIndex};
           [~, HTest, error, Pxx, f, t] = nmf_eval_v2(mixtureSignal(:,:,nbitsIndex,loopIndex), paramsNMF2);
            
            for JNRIndex = 1:length(paramsNMF1.JNRVector)
                for i = 1:2
                    S = (paramsNMF2.W0(:,(i-1)*paramsNMF2.numberOfSources/2 +1:(i*paramsNMF2.numberOfSources/2),nbitsIndex) * ...
                        HTest{1,JNRIndex}((i-1)*paramsNMF2.numberOfSources/2 +1:(i*paramsNMF2.numberOfSources/2),:) ./ (paramsNMF2.W0(:,:,nbitsIndex)*HTest{1,JNRIndex})).*Pxx{1,JNRIndex};
                    
                    V = paramsNMF2.W0(:,(i-1)*paramsNMF2.numberOfSources/2 +1:(i*paramsNMF2.numberOfSources/2),nbitsIndex) * ...
                        HTest{1,JNRIndex}((i-1)*paramsNMF2.numberOfSources/2 +1:(i*paramsNMF2.numberOfSources/2),:) .* exp(1j*(angle(PxxGPS)+1*randn(size(PxxGPS))));
                    xHat(:,i,JNRIndex,nbitsIndex, bandwidthIndex, loopIndex) = istft(S, paramsNMF1.fs, 'Window', paramsNMF1.window, 'OverlapLength', paramsNMF1.overlap, 'FFTLength', paramsNMF1.nfft);
                end
            end
        end
    end
end

v = angle(PxxGPS);

save(['.' filesep 'data' filesep 'nmf_testing_33_1.mat'], 'xHat', 'nbits', 'JNRVector', 'v', '-v7.3');

rmpath(['.' filesep 'data']);
rmpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);