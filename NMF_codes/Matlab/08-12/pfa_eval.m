clear;
clc;
close all;

addpath(['..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);

load sim_params_1.mat;

random_state = 42;

params.fs = paramsSignal.Freqsamp;
params.nfft = 128;
params.nperseg = 128;
params.overlap = params.nperseg-1;
params.hop_size = params.nperseg - params.overlap;
params.numberOfSources = 1;
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 10000;
params.tolChange = 1e-6;
params.tolError = 1e-6;
params.repetitions = 1;
params.JNRVector = [-17];

bandwidthVector = 10.72e6;
periodVector = 8.72e-6;

rng(random_state)

numberOfRawSamples = 4096;
silenceSamples = round(20e-6*params.fs);
thresholdVector = 0.1:0.05:0.9;
window_median_length_vector = 51:50:401;
monteCarloLoops = 100;
noisePowerVector = [0.7209 2.2798 7.2092 22.7975 72.09];
outputLength = (numberOfRawSamples + silenceSamples*2 - params.nperseg + 1)/(params.nperseg - params.overlap);
detection_res = zeros(monteCarloLoops, length(noisePowerVector), length(thresholdVector), length(window_median_length_vector), outputLength);

for loopIndex = 1:monteCarloLoops
    loopIndex
    noise = randn(numberOfRawSamples + silenceSamples*2, 1) + 1j*randn(numberOfRawSamples + silenceSamples*2, 1);
    noisePower = pow_eval(noise);
    noise2 = zeros(length(noise), length(noisePowerVector));
    
    for noisePowerIndex = 1:length(noisePowerVector)
        noise2(:,noisePowerIndex) = noise*sqrt(db2pow(noisePowerVector(noisePowerIndex))/noisePower);
    end
    
    [W, ~, ~, PxxAux, ~, ~] = nmf_eval_v2(noise2, params);
    
    for JNRIndex = 1:length(noisePowerVector)
        
        inputNMF = abs(PxxAux{1, JNRIndex}).^2;
        
        inputNMF = inputNMF - mean(inputNMF);
        inputNMF = inputNMF.*sqrt(1./var(inputNMF));
        inputNMFAux = sqrt(sum(inputNMF.*inputNMF)) + eps;
        inputNMFNormalised = inputNMF./inputNMFAux;
        
        WNormalised = W{1, JNRIndex}(:,1) - mean(W{1, JNRIndex}(:,1));
        WNormalised = WNormalised.*sqrt(1./var(WNormalised));
        WNormalised = WNormalised ./ (norm(WNormalised) + eps);
        
        output = inputNMFNormalised.'*WNormalised;
        for thresholdIndex = 1:length(thresholdVector)
            for window_median_length_index = 1:length(window_median_length_vector)
                detection_res(loopIndex, JNRIndex, thresholdIndex, window_median_length_index, :) = ...
                    detection_eval(output, thresholdVector(thresholdIndex), window_median_length_vector(window_median_length_index));
            end
        end
    end
end

save(['..' filesep '.' filesep 'data' filesep '08-12' filesep 'results05.mat'], 'detection_res', '-v7.3');

rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep '.' filesep 'Sigtools' filesep])