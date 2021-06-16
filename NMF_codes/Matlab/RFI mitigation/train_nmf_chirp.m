clear;
clc;
close all;
    
addpath(['..' filesep 'Sigtools' filesep])
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])

load(['..' filesep  'signalsGeneration' filesep 'sim_params' filesep 'sim_params_3.mat']);

monteCarloLoops = 1;
SNR = -20;

params.JNRVector = 0;
params.fs = paramsSignal.Freqsamp;
params.nfft = 256;
params.nperseg = 256;
params.overlap = params.nperseg - 1;
params.hop_size = params.nperseg - params.overlap;
params.window = kaiser(params.nperseg, 60);
params.specType = 'power';
params.numberOfSources = 60;
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 500;
params.tolChange = 1e-6;
params.tolError = 1e-6;
params.repetitions = 1;
params.transpose = false;
params.verbose = true;
paramsSignal.FreqDopp = 1.5e3;

bandwidthVector = (2:6:14)*1e6;
periodVector = 8.62e-6;
paramsSignal.Noneperiod = round(periodVector*params.fs);                   % number of samples with a sweep time

numberOfRawSamples = floor(paramsSignal.Freqsamp*paramsSignal.Intetime);
totalSamples = floor(paramsSignal.Freqsamp*paramsSignal.Intetime);

GPSSignals = GPSGen(paramsSignal);
GPSSignals = GPSSignals(1:numberOfRawSamples,:);
GPSSignalsPower = pow_eval(GPSSignals);

W0 = cell(1, length(bandwidthVector));

for loopIndex = 1:monteCarloLoops
    loopIndex
    noise = randn(totalSamples, 1) + 1j*randn(totalSamples, 1);
    noisePower = pow_eval(noise);

    GPSSignalsAux = GPSSignals;
    GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
    mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
    [WSignal, HSignal, errorSignalTrain, PxxSignal, f, t] = nmf_eval_v2(mixtureGPS, params);
    
    for bandwidthIndex = 1:length(bandwidthVector)
        
        paramsSignal.IFmin = -bandwidthVector(bandwidthIndex)/2;                                                  % start frequency
        paramsSignal.IFmax = bandwidthVector(bandwidthIndex)/2;                    % end frequency
        paramsSignal.foneperiod(1:paramsSignal.Noneperiod) = linspace(paramsSignal.IFmin, paramsSignal.IFmax, paramsSignal.Noneperiod);
        
        [interferenceSignal, ~] = interferenceGen(paramsSignal);
        interferenceSignal = interferenceSignal(1:numberOfRawSamples);
        interferenceSignalPower = pow_eval(interferenceSignal);
      
        [WInterf, HInterf, errorInterfTrain, PxxInterf, ~, ~] = nmf_eval_v2(interferenceSignal, params);
        W0{1,bandwidthIndex} = [WInterf{1,1} WSignal{1,1}];
        
    end
     
end

save(['.' filesep 'data' filesep 'nmf_training_3.mat'], 'W0', 'params');

rmpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep]);