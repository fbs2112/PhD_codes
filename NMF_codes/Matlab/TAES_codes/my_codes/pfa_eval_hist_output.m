clear;
clc;
close all;

addpath(['..' filesep '..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep '..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['..' filesep '..' filesep  'signalsGeneration' filesep]);
addpath(['..' filesep '..' filesep  'signalsGeneration' filesep 'sim_params']);

load sim_params_1.mat;

random_state = 42;
rng(random_state)

params.fs = paramsSignal.Freqsamp;
params.nfft = 64;
params.nperseg = 64;
params.overlap = params.nperseg - 1;

params.hop_size = params.nperseg - params.overlap;
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
monteCarloLoops = 1000;

GPSSignals = GPSGen(paramsSignal);
GPSSignals = GPSSignals(1:numberOfRawSamples,:);
GPSSignalsPower = pow_eval(GPSSignals);
outputLength = floor((numberOfRawSamples - params.overlap)/(params.hop_size));

output_res = zeros(monteCarloLoops, 1);

for loopIndex = 1:monteCarloLoops
    loopIndex
    noise = randn(totalSamples, 1) + 1j*randn(totalSamples, 1);
    noisePower = pow_eval(noise);
   
    GPSSignalsAux = GPSSignals;
    GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
    mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
    mixtureSignal = mixtureGPS ;
    
    [W, ~, ~, PxxAux, ~, ~] = nmf_eval_v2(mixtureSignal, params);
    
    inputNMF = abs(PxxAux{1, 1}).^2;
    
    inputNMF = inputNMF - mean(inputNMF);
    inputNMF = inputNMF.*sqrt(1./var(inputNMF));
    inputNMFAux = sqrt(sum(inputNMF.*inputNMF)) + eps;
    inputNMFNormalised = inputNMF./inputNMFAux;
    
    WNormalised = W{1, 1}(:,1) - mean(W{1, 1}(:,1));
    WNormalised = WNormalised.*sqrt(1./var(WNormalised));
    WNormalised = WNormalised ./ (norm(WNormalised) + eps);
    
    output = inputNMFNormalised.'*WNormalised;
    output_res(loopIndex) = median(output);
   
end

save(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results' filesep 'results46.mat'], 'output_res', '-v7.3');

rmpath(['..' filesep '..' filesep '.' filesep 'Sigtools' filesep])
rmpath(['..' filesep '..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep '..' filesep  'signalsGeneration' filesep]);
rmpath(['..' filesep '..' filesep  'signalsGeneration' filesep 'sim_params']);