clear;
clc;
close all;

addpath(['..' filesep 'Sigtools' filesep])
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
addpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['.' filesep 'data']);

load nmf_training_19.mat;
load sim_params_3.mat;

monteCarloLoops = 100;
SNR = -25;
nbits = 0;

params.JNRVector = [-5 0 10 30 50];
% params.JNRVector = [30 50];

JNRVector = params.JNRVector; 
params.fs = paramsSignal.Freqsamp;
params.nfft = 256;
params.nperseg = 256;
params.overlap = params.nperseg - 1;
params.hop_size = params.nperseg - params.overlap;
params.window = ones(params.nperseg, 1);
params.specType = 'power';
params.numberOfSources = 10;
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 500;
params.tolChange = 1e-3;
params.tolError = 1e-3;
params.repetitions = 1;
params.verbose = false;

numberOfRawSamples = floor(paramsSignal.Freqsamp*paramsSignal.Intetime);
delay = 10e-6;
totalSamples = numberOfRawSamples;

initialFrequency = 2e6;
bandwidthVector = 2e6;
periodVector = 8.62e-6;

paramsSignal.Noneperiod = round(periodVector*params.fs);                   % number of samples with a sweep time
paramsSignal.IFmin = initialFrequency;                                                  % start frequency
paramsSignal.IFmax = bandwidthVector + initialFrequency;                    % end frequency
paramsSignal.foneperiod(1:paramsSignal.Noneperiod) = linspace(paramsSignal.IFmin, paramsSignal.IFmax, paramsSignal.Noneperiod);
paramsSignal.Initphase = 0;

interferenceSignal = interferenceGen(paramsSignal);
interferenceSignal = interferenceSignal(1:numberOfRawSamples);

Timeofthisloop = 0:totalSamples-1;
Carrphase = mod(2*pi*(paramsSignal.FreqDopp)*Timeofthisloop/paramsSignal.Freqsamp,2*pi);
Carrier = exp(1i*Carrphase).';

interferenceSignal = interferenceSignal.*Carrier;
interferenceSignalPower = pow_eval(interferenceSignal);
  
mixtureSignal = zeros(totalSamples, length(params.JNRVector), length(nbits), monteCarloLoops);
xHat = zeros(totalSamples, 2, length(params.JNRVector), length(nbits), monteCarloLoops);
xHatSemi = zeros(totalSamples, 2, length(params.JNRVector), length(nbits), monteCarloLoops);

paramsSignal.FreqDopp = 1.2e3;
paramsSignal.numberOfGPSSignals = 1;

GPSSignals = GPSGen(paramsSignal);
GPSSignals = GPSSignals(1:numberOfRawSamples,:);
GPSSignals = [GPSSignals(end - round(delay*params.fs)+1:end,:);GPSSignals(1:end - round(delay*params.fs),:)]; % Introducing artificial code delay
GPSSignalsPower = pow_eval(GPSSignals);

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

    for nbitsIndex = 1:length(nbits)
        nbitsIndex
        
        for JNRIndex = 1:length(params.JNRVector)
            GPSSignalsAux = GPSSignals;
            interferenceSignalAux = interferenceSignal;
            GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
            mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
            interferenceSignalAux = interferenceSignalAux*sqrt(noisePower*10^(params.JNRVector(JNRIndex)/10)/interferenceSignalPower);
            mixtureSignal(:,JNRIndex,nbitsIndex,loopIndex) = mixtureGPS + interferenceSignalAux;
%             mixtureSignal(:,JNRIndex,nbitsIndex,loopIndex) = quantise_gps(mixtureSignal(:,JNRIndex,nbitsIndex,loopIndex), nbits(nbitsIndex), noiseVar);
        end
        params.init = 'random';
        params.transform = true;
        params.semi = false;

        [WTestAux, ~, ~, ~, ~, ~] = nmf_eval_v2(mixtureSignal(:,:,nbitsIndex,loopIndex), params);
        
        params.numberOfSources = params.numberOfSources*2;
        params.init = 'custom';
        params.transform = false;
        params.semi = true;
        %--------Semi supervised NMF
        for idx = 1:length(params.JNRVector)
            WOSemi(:,(idx-1) * params.numberOfSources + 1:idx * params.numberOfSources) = [WTestAux{1,idx} W0(:,params.numberOfSources/2+1:end,nbitsIndex)];
        end
        params.W0 = WOSemi;
        [~, HTest, error, Pxx, f, t] = nmf_eval_v2(mixtureSignal(:,:,nbitsIndex,loopIndex), params);
        
        for JNRIndex = 1:length(params.JNRVector)
            wAux = WOSemi(:,(JNRIndex-1) * params.numberOfSources + 1:JNRIndex * params.numberOfSources);
            for i = 1:2
                S = (wAux(:,(i-1)*params.numberOfSources/2 +1:(i*params.numberOfSources/2),nbitsIndex) * ...
                    HTest{1,JNRIndex}((i-1)*params.numberOfSources/2 +1:(i*params.numberOfSources/2),:) ./ (wAux(:,:,nbitsIndex)*HTest{1,JNRIndex})).*Pxx{1,JNRIndex};
                
                xHat(:,i,JNRIndex, nbitsIndex,loopIndex) = istft(S, params.fs, 'Window', params.window, 'OverlapLength', params.overlap, 'FFTLength', params.nfft);
            end
        end
        params.numberOfSources = params.numberOfSources/2;

    end
end

save(['.' filesep 'data' filesep 'nmf_testing_31.mat'], 'xHat', 'mixtureSignal', 'nbits', 'JNRVector');

rmpath(['.' filesep 'data']);
rmpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);