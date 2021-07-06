clear;
clc;
close all;

addpath(['..' filesep 'Sigtools' filesep])
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])

load(['.' filesep 'data' filesep 'nmf_training_5.mat']);
load(['..' filesep  'signalsGeneration' filesep 'sim_params' filesep 'sim_params_3.mat']);

monteCarloLoops = 100;
SNR = -20;
numberOfSources = 2;

params.JNRVector = 0:5:30;

JNRVector = params.JNRVector;
paramsNMF = params;
paramsNMF.numberOfSources = paramsNMF.numberOfSources*numberOfSources;
paramsNMF.init = 'custom';
paramsNMF.transform = false;
paramsNMF.verbose = false;

numberOfZerosVector = params.fs*(1e-3)*[0];

numberOfRawSamples = floor(paramsSignal.Freqsamp*paramsSignal.Intetime);
delay = 10e-6;
totalSamples = numberOfRawSamples;

bandwidthVector = [2 8 14]*1e6;
periodVector = 8.62e-6;

paramsSignal.Noneperiod = round(periodVector*paramsNMF.fs);                   % number of samples with a sweep time
paramsSignal.Initphase = 0;
paramsSignal.FreqDopp = 1e3;

GPSSignals = GPSGen(paramsSignal);
GPSSignals = GPSSignals(1:numberOfRawSamples,:);
GPSSignals = [GPSSignals(end - round(delay*paramsNMF.fs)+1:end,:);GPSSignals(1:end - round(delay*paramsNMF.fs),:)]; % Introducing artificial code delay
GPSSignalsPower = pow_eval(GPSSignals);
edgeZeros = 100;

mixtureSignal = zeros(totalSamples + 2*edgeZeros, length(paramsNMF.JNRVector), monteCarloLoops);
xHat = zeros(totalSamples, numberOfSources, length(paramsNMF.JNRVector), length(numberOfZerosVector), length(bandwidthVector), monteCarloLoops);

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
    
    for bandwidthIndex = 1:length(bandwidthVector)
        bandwidthIndex
        paramsSignal.IFmin = -bandwidthVector(bandwidthIndex)/2;                                                  % start frequency
        paramsSignal.IFmax = bandwidthVector(bandwidthIndex)/2;                    % end frequency
        paramsSignal.foneperiod(1:paramsSignal.Noneperiod) = linspace(paramsSignal.IFmin, paramsSignal.IFmax, paramsSignal.Noneperiod);
        
        [interferenceSignal, ~] = interferenceGen(paramsSignal);
        interferenceSignal = interferenceSignal(1:numberOfRawSamples);
        interferenceSignalPower = pow_eval(interferenceSignal);
        
        for numberOfZerosIndex = 1:length(numberOfZerosVector)
            for JNRIndex = 1:length(paramsNMF.JNRVector)
                GPSSignalsAux = GPSSignals;
                interferenceSignalAux = interferenceSignal;
                interferenceSignalAux(1:numberOfZerosVector(numberOfZerosIndex)) = 0;
                
                GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
                mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
                interferenceSignalAux = interferenceSignalAux*sqrt(noisePower*10^(paramsNMF.JNRVector(JNRIndex)/10)/interferenceSignalPower);
                mixtureSignalZeros = mixtureGPS + interferenceSignalAux;
                mixtureSignalZeros = [zeros(edgeZeros, 1); mixtureSignalZeros; zeros(edgeZeros, 1)];
                mixtureSignal(:,JNRIndex,loopIndex) = mixtureSignalZeros;
            end
            paramsNMF.W0 = W0{1,bandwidthIndex}(:,1:params.numberOfSources*numberOfSources);
            
            [W, H, error, Pxx, f, t] = nmf_eval_v2(mixtureSignal(:,:,loopIndex), paramsNMF);
            
            for JNRIndex = 1:length(paramsNMF.JNRVector)
                for i = 1:numberOfSources
                    S = (W{1,JNRIndex}(:,(i-1)*paramsNMF.numberOfSources/numberOfSources +1:(i*paramsNMF.numberOfSources/numberOfSources)) * ...
                        H{1,JNRIndex}((i-1)*paramsNMF.numberOfSources/numberOfSources +1:(i*paramsNMF.numberOfSources/numberOfSources),:) ./ (W{1,JNRIndex}*H{1,JNRIndex})).*Pxx{1,JNRIndex};
                    
                    if paramsNMF.transpose
                        S = S.';
                    end
                    
                    if ~isfield(params, 'tf')
                        xHatAux = istft(S, paramsNMF.fs, 'Window', paramsNMF.window, 'OverlapLength', paramsNMF.overlap, 'FFTLength', paramsNMF.nfft);
                    elseif isfield(params, 'tf') && strcmp(params.tf, 'fsst')
                        xHatAux = ifsst(S, params.window);
                    end
                    xHat(:,i,JNRIndex,numberOfZerosIndex,bandwidthIndex,loopIndex) = xHatAux(edgeZeros+1:end-edgeZeros); %removing zeros at the edge
                    
                end
            end
        end
        
    end
end

xHat = single(xHat);

save(['.' filesep 'data' filesep 'nmf_testing_12.mat'], 'xHat', 'JNRVector', '-v7.3');

rmpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep]);