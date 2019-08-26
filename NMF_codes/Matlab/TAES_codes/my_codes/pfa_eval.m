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
params.nfft = 256;
params.nperseg = 128;
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
thresholdVector = 0.1:0.05:0.9;
window_median_length_vector = 51:50:401;
monteCarloLoops = 1000;

GPSSignals = GPSGen(paramsSignal);
GPSSignals = GPSSignals(1:numberOfRawSamples,:);
GPSSignalsPower = pow_eval(GPSSignals);
    
outputLength = floor((totalSamples - params.overlap)/(params.nperseg - params.overlap));
detection_res = zeros(monteCarloLoops, length(thresholdVector), length(window_median_length_vector), outputLength);


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
    
    for thresholdIndex = 1:length(thresholdVector)
        for window_median_length_index = 1:length(window_median_length_vector)
            detection_res(loopIndex, thresholdIndex, window_median_length_index, :) = detection_eval(output, thresholdVector(thresholdIndex), ...
                window_median_length_vector(window_median_length_index));
        end
    end
end

if isunix
    save(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep 'pfa_results' filesep 'results07.mat'], 'detection_res', '-v7.3');
else
end

rmpath(['..' filesep '..' filesep '.' filesep 'Sigtools' filesep])
rmpath(['..' filesep '..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep '..' filesep  'signalsGeneration' filesep]);
rmpath(['..' filesep '..' filesep  'signalsGeneration' filesep 'sim_params']);