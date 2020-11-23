clear;
clc;
close all;

frequency = '1152';
% fileName = '2018-09-01-09_52_55_0000020535312384.000000'; %train
fileName = '2018-09-01-09_52_55_0000000000000000.000000'; %test

load(['.' filesep 'data' filesep 'template_ADSB.mat']);
load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep fileName '_labels.mat']);
load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep fileName '_labels_train_test.mat']);

params.fs = 128e6;
params.nfft = 256;
params.nperseg = 256;
params.overlap = params.nperseg/2;
params.hop_size = params.nperseg - params.overlap;
params.window = ones(params.nperseg, 1);
params.JNRVector = 0;
thresholdVector = 0:0.0005:1;
window_median_length_vector = 0;

% falseLabels = find(interferenceDetFlag == 0 | interferenceDetFlag == 1 | interferenceDetFlag == 4);
% monteCarloLoops = length(falseLabels);
monteCarloLoops = length(falseLabelsTrain);

detection_res = zeros(monteCarloLoops, length(thresholdVector), length(window_median_length_vector));
feature = zeros(monteCarloLoops, 1);
for loopIndex = 1:monteCarloLoops
    loopIndex
    load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep fileName '_' num2str(falseLabelsTrain(loopIndex)) '.mat']);  
    [Pxx, f, t] = spectrogram(parkesSignalA, params.window, params.overlap, params.nfft, params.fs);
    PxxAux{1,1} = Pxx;
    for JNRIndex = 1:length(params.JNRVector)
        
        inputNMF = abs(PxxAux{1, JNRIndex}).^2;
        inputNMF = inputNMF - mean(inputNMF);
        inputNMF = inputNMF.*sqrt(1./var(inputNMF));
        inputNMFAux = sqrt(sum(inputNMF.*inputNMF)) + eps;
        inputNMFNormalised = inputNMF./inputNMFAux;
        
        templateNormalised = template - mean(template);
        templateNormalised = templateNormalised.*sqrt(1./var(templateNormalised));
        templateNormalised = templateNormalised ./ (norm(templateNormalised) + eps);
        
        output = max((inputNMFNormalised.'*templateNormalised));
        feature(loopIndex) = output;
        for thresholdIndex = 1:length(thresholdVector)
            for window_median_length_index = 1:length(window_median_length_vector)
                detection_res(loopIndex, thresholdIndex, window_median_length_index) = ((output> thresholdVector(thresholdIndex)));
            end
        end
    end    
end

save(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_57.mat'], 'detection_res', 'thresholdVector', 'feature', '-v7.3');