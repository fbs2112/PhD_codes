clear;
clc;
close all;


matNumber = 53;
frequency = '1152';
load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_' num2str(matNumber) '.mat']);

% thresholdVector = 0:0.005:1;
% thresholdVectorAux = 0:0.005:1;
% thresholdVector = linspace(380, 420, length(thresholdVectorAux));
window_median_length_vector = 0;
fp = zeros(size(detection_res,1), length(thresholdVector));
tn = zeros(size(detection_res,1), length(thresholdVector));

for k = 1:length(thresholdVector)
    for j = 1:length(window_median_length_vector)
        x = squeeze(detection_res(:,k));
        fp(:,k) = x;
        tn(:,k) = 1 - fp(:,k) ;
    end
end

fpr = fp./(fp+tn);
averageFpr = squeeze(mean(fpr, 1));
stdFpr = squeeze(std(fpr, [], 1));

figure;
plot(thresholdVector, averageFpr)

ylabel('Probability of false alarm');
xlabel('$\bar{\gamma}$');
xlim([min(thresholdVector) max(thresholdVector)]);
grid on;

figure;
semilogy(thresholdVector, averageFpr)

ylabel('Probability of false alarm');   
xlabel('$\bar{\gamma}$');
xlim([min(thresholdVector) max(thresholdVector)]);
grid on;

save(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'pfa_det_' num2str(matNumber) '.mat'], 'averageFpr', 'stdFpr', '-v7.3');
