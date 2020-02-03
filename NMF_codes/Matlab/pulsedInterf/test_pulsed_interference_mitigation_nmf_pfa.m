clear;
clc;
close all;

addpath(['..' filesep 'Sigtools' filesep])
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
addpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['.' filesep 'data']);

load nmf_training_04.mat;
load sim_params_2.mat;

monteCarloLoops = 100;
SNR = -25;
alpha = 4.5e11;
deltaT = 12e-6;

params.fs = paramsSignal.Freqsamp;
params.nfft = 256;
params.nperseg = 256;
params.overlap = params.nperseg - 1;
params.hop_size = params.nperseg - params.overlap;
params.window = ones(params.nperseg, 1);
params.specType = 'power';
params.numberOfSources = 5;
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 1000;
params.tolChange = 1e-6;
params.tolError = 1e-6;
params.repetitions = 1;
params.JNRVector = 0;

numberOfRawSamples = floor(paramsSignal.Freqsamp*paramsSignal.Intetime);
totalSamples = numberOfRawSamples;

paramsSignal.numberOfGPSSignals = 1;
paramsSignal.FreqDopp = 3e3;
params.numberOfSources = params.numberOfSources*2;
params.init = 'custom';
params.W0 = W0;
params.transform = false;

mixtureSignal = zeros(totalSamples, monteCarloLoops);
xHat = zeros(length(mixtureSignal), 2, 1, monteCarloLoops);

for loopIndex = 1:monteCarloLoops
    loopIndex
    noise = randn(totalSamples, 1) + 1j*randn(totalSamples, 1);
    
    mixtureSignal(:,loopIndex) = noise;
    [WTest, HTest, errorTest, PxxTest, f, t] = nmf_eval_v2(mixtureSignal(:,loopIndex), params);
    S = zeros([size(PxxTest{1,1}) 2]);
    
    for i = 1:2
        S(:,:,i) = (W0(:,(i*params.numberOfSources/2) - (params.numberOfSources/2 -1):(i*params.numberOfSources/2)) * ...
            HTest{1,1}((i*params.numberOfSources/2) - (params.numberOfSources/2 -1):(i*params.numberOfSources/2),:)./ (W0*HTest{1,1})).*PxxTest{1,1};
        
        xHat(:,i,1,loopIndex) = istft(S(:,:,i), params.fs, 'Window', params.window, 'OverlapLength', params.overlap, 'FFTLength', params.nfft);
        
    end
end

save(['.' filesep 'data' filesep 'nmf_testing_12.mat'], 'xHat', 'mixtureSignal');

rmpath(['.' filesep 'data']);
rmpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);