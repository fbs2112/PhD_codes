clear;
clc;
close all;

addpath(['..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['..' filesep  'signalsGeneration' filesep]);
addpath(['..' filesep  'signalsGeneration' filesep 'sim_params']);

load sim_params_1.mat;

warning('off','all')

params.fs = paramsSignal.Freqsamp;
params.nfft = 32;
params.nperseg = 32;
params.overlap = 0;
params.hop_size = params.nperseg - params.overlap;
params.window = ones(params.nperseg, 1);
params.specType = 'power';
params.numberOfSources = 1;
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 10000;
params.tolChange = 1e-6;
params.tolError = 1e-6;
params.repetitions = 1;
params.JNRVector = 0;
SNR = -25;

numberOfRawSamples = 4096;
totalSamples = numberOfRawSamples;
thresholdVector = 0:0.05:0.5;
monteCarloLoops = 100;

timeBins = floor((totalSamples - params.overlap)/(params.nperseg - params.overlap));
detection_res = zeros(monteCarloLoops, timeBins, length(thresholdVector));
pvalue = zeros(monteCarloLoops, timeBins);
kvalue = zeros(monteCarloLoops, timeBins);

for loopIndex = 1:monteCarloLoops
    loopIndex
    noise = randn(totalSamples, 1) + 1j*randn(totalSamples, 1);
    noisePower = pow_eval(noise);
    
    GPSSignals = GPSGen(paramsSignal);
    GPSSignals = GPSSignals(1:numberOfRawSamples,:);
    GPSSignalsPower = pow_eval(GPSSignals);
    
    GPSSignalsAux = GPSSignals;
    GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
    mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
    mixtureSignal = mixtureGPS;
    
    [detection_res(loopIndex, :, :), pvalue(loopIndex, :), ~] = ...
        deteTimeSlice(mixtureSignal, params, thresholdVector);
end

save(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'statistical_data' filesep ...
        filesep 'pfa_results' filesep 'results1.mat'], 'detection_res', '-v7.3');

rmpath(['..' filesep '.' filesep 'Sigtools' filesep])
rmpath(['..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep  'signalsGeneration' filesep]);
rmpath(['..' filesep  'signalsGeneration' filesep 'sim_params']);

warning('on','all')