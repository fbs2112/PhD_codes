clear;
clc;
close all;

addpath(['..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['.' filesep 'data']);

load ADSB_label.mat;

params.fs = 128e6;
params.nfft = 256;
params.nperseg = 256;
params.overlap = params.nperseg/2;
params.hop_size = params.nperseg - params.overlap;
params.window = ones(params.nperseg, 1);
params.specType = 'power';
params.numberOfSources = 2;
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 10000;
params.tolChange = 1e-6;
params.tolError = 1e-6;
params.repetitions = 1;
params.JNRVector = 0;

thresholdVector = 0:0.005:1;
window_median_length_vector = 0;

falseLabels = find(~interferenceDetFlag);
monteCarloLoops = length(falseLabels);

detection_res = zeros(monteCarloLoops, length(thresholdVector), length(window_median_length_vector));
for loopIndex = 1:monteCarloLoops
    loopIndex
    load(['.' filesep 'data' filesep 'dataParkes_' num2str(falseLabels(loopIndex)) '.mat']);   
    [W, ~, ~, PxxAux, f, t] = nmf_eval_detection(parkesSignal, params);
    
    for JNRIndex = 1:length(params.JNRVector)
        
        inputNMF = abs(PxxAux{1, JNRIndex}).^2;
        inputNMF = inputNMF - mean(inputNMF);
        inputNMF = inputNMF.*sqrt(1./var(inputNMF));
        inputNMFAux = sqrt(sum(inputNMF.*inputNMF)) + eps;
        inputNMFNormalised = inputNMF./inputNMFAux;
        
        WNormalised1 = W{1, JNRIndex}(:,1) - mean(W{1, JNRIndex}(:,1));
        WNormalised1 = WNormalised1.*sqrt(1./var(WNormalised1));
        WNormalised1 = WNormalised1 ./ (norm(WNormalised1) + eps);
        
        WNormalised2 = W{1, JNRIndex}(:,2) - mean(W{1, JNRIndex}(:,2));
        WNormalised2 = WNormalised2.*sqrt(1./var(WNormalised2));
        WNormalised2 = WNormalised2 ./ (norm(WNormalised2) + eps);
        
        output1 = ((inputNMFNormalised.'*WNormalised1));
        output2 = ((inputNMFNormalised.'*WNormalised2));
        for thresholdIndex = 1:length(thresholdVector)
            for window_median_length_index = 1:length(window_median_length_vector)
                detection_res(loopIndex, thresholdIndex, window_median_length_index) = ...
                    (detection_eval_2(output, thresholdVector(thresholdIndex)));
            end
        end
    end    
end

% save(['.' filesep 'data' filesep 'results_det_2.mat'], 'detection_res', '-v7.3');

rmpath(['.' filesep 'data']);
rmpath(['..' filesep '.' filesep 'Sigtools' filesep])
rmpath(['..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])