clear;
clc;
close all;

frequency = '1152';
% fileName = '2018-09-01-09_52_55_0000020535312384.000000'; %train
fileName = '2018-09-01-09_52_55_0000000000000000.000000'; %test

load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep fileName '_labels.mat']);
load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep fileName '_labels_train_test.mat']);

params.fs = 128e6;
params.nfft = 256;
params.nperseg = 256;
params.overlap = params.nperseg/2;
params.hop_size = params.nperseg - params.overlap;
params.window = ones(params.nperseg, 1);
params.JNRVector = 0;

thresholdVectorAux = 0:0.0005:1;
thresholdVector = 200:0.5:3000;

window_median_length_vector = 0;

% trueLabels = find(interferenceDetFlag == 2 | interferenceDetFlag == 3);
monteCarloLoops = length(trueLabelsTest);

detection_res = zeros(monteCarloLoops, length(thresholdVector), length(window_median_length_vector));
feature = zeros(monteCarloLoops, 1);

for loopIndex = 1:monteCarloLoops
    loopIndex
    
    load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep fileName '_' num2str(trueLabelsTest(loopIndex)) '.mat']);  
    parkesSignalBuffer = buffer(parkesSignalA, params.nfft, params.overlap, 'nodelay');
    
    detAux = max(abs(goertzel(parkesSignalBuffer, 132, 1)));
    feature(loopIndex) = detAux;
    for thresholdIndex = 1:length(thresholdVector)
        for window_median_length_index = 1:length(window_median_length_vector)
            detection_res(loopIndex, thresholdIndex, window_median_length_index) = (detAux > thresholdVector(thresholdIndex));
        end
    end
    
end

save(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_54.mat'], 'detection_res', 'thresholdVector', 'feature', '-v7.3');
