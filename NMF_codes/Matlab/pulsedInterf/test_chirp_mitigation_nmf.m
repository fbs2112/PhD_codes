clear;
clc;
close all;

addpath(['..' filesep 'Sigtools' filesep])
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
addpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['.' filesep 'data']);

load nmf_training_30.mat;
load sim_params_3.mat;

monteCarloLoops = 100;
SNR = -20;
nbits = 0;
numberOfSources = 2;

params.JNRVector = 0:5:30;

JNRVector = params.JNRVector;
paramsNMF = params;
paramsNMF.numberOfSources = paramsNMF.numberOfSources*numberOfSources;
paramsNMF.init = 'custom';
paramsNMF.transform = false;

numberOfZerosVector = params.fs*(1e-3)*[0 0.25 0.5];

numberOfRawSamples = floor(paramsSignal.Freqsamp*paramsSignal.Intetime);
delay = 10e-6;
totalSamples = numberOfRawSamples;

bandwidthVector = (2:6:14)*1e6;
periodVector = 8.62e-6;

paramsSignal.Noneperiod = round(periodVector*paramsNMF.fs);                   % number of samples with a sweep time
paramsSignal.Initphase = 0;
paramsSignal.FreqDopp = 1e3;

GPSSignals = GPSGen(paramsSignal);
GPSSignals = GPSSignals(1:numberOfRawSamples,:);
GPSSignals = [GPSSignals(end - round(delay*paramsNMF.fs)+1:end,:);GPSSignals(1:end - round(delay*paramsNMF.fs),:)]; % Introducing artificial code delay
GPSSignalsPower = pow_eval(GPSSignals);

mixtureSignal = zeros(totalSamples, length(paramsNMF.JNRVector), length(nbits), monteCarloLoops);
xHat = zeros(totalSamples, numberOfSources, length(paramsNMF.JNRVector), length(numberOfZerosVector), length(nbits), length(bandwidthVector), monteCarloLoops);

for loopIndex = 1:monteCarloLoops
    loopIndex
    
    if paramsSignal.FreqDopp ~= 0
        noise = randn(totalSamples, 1) + 1j*randn(totalSamples, 1)/20;
        noiseVar = 2;
    else
        noise = randn(totalSamples, 1);
        noiseVar = 1;
    end
    noisePower = pow_eval(noise);
    
    for bandwidthIndex = 1:length(bandwidthVector)
        bandwidthIndex
        paramsSignal.IFmin = -bandwidthVector(bandwidthIndex)/2;                                                  % start frequency
        paramsSignal.IFmax = bandwidthVector(bandwidthIndex)/2;                    % end frequency
        paramsSignal.foneperiod(1:paramsSignal.Noneperiod) = linspace(paramsSignal.IFmin, paramsSignal.IFmax, paramsSignal.Noneperiod);
        
        [interferenceSignal, ~] = interferenceGen(paramsSignal);
        interferenceSignal = interferenceSignal(1:numberOfRawSamples);
        interferenceSignalPower = pow_eval(interferenceSignal);
        
        for nbitsIndex = 1:length(nbits)
            for numberOfZerosIndex = 1:length(numberOfZerosVector)
                for JNRIndex = 1:length(paramsNMF.JNRVector)
                    GPSSignalsAux = GPSSignals;
                    interferenceSignalAux = interferenceSignal;
                    interferenceSignalAux(1:numberOfZerosVector(numberOfZerosIndex)) = 0;
                    
                    GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
                    mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
                    interferenceSignalAux = interferenceSignalAux*sqrt(noisePower*10^(paramsNMF.JNRVector(JNRIndex)/10)/interferenceSignalPower);
                    mixtureSignal(:,JNRIndex,nbitsIndex,loopIndex) = mixtureGPS + interferenceSignalAux;
                end
                paramsNMF.W0 = W0{1,bandwidthIndex}(:,1:params.numberOfSources*numberOfSources);
                                
                [W, H, error, Pxx, f, t] = nmf_eval_v2(mixtureSignal(:,:,nbitsIndex,loopIndex), paramsNMF);
                
                for JNRIndex = 1:length(paramsNMF.JNRVector)
                    for i = 1:numberOfSources
                        S = (W{1,JNRIndex}(:,(i-1)*paramsNMF.numberOfSources/numberOfSources +1:(i*paramsNMF.numberOfSources/numberOfSources),nbitsIndex) * ...
                            H{1,JNRIndex}((i-1)*paramsNMF.numberOfSources/numberOfSources +1:(i*paramsNMF.numberOfSources/numberOfSources),:) ./ (W{1,JNRIndex}*H{1,JNRIndex})).*Pxx{1,JNRIndex};
                        
                        if paramsNMF.transpose
                            S = S.';
                        end
                        
                        xHat(:,i,JNRIndex,numberOfZerosIndex,nbitsIndex,bandwidthIndex,loopIndex) = istft(S, paramsNMF.fs, 'Window', paramsNMF.window, 'OverlapLength', paramsNMF.overlap, 'FFTLength', paramsNMF.nfft);
                        
                    end
                end
            end
        end
    end
end

xHat = single(xHat);

save(['.' filesep 'data' filesep 'nmf_testing_35.mat'], 'xHat', 'nbits', 'JNRVector', '-v7.3');

rmpath(['.' filesep 'data']);
rmpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);