clear;
clc;
% close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '..' filesep '.' filesep 'Misc'])
addpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results']);  

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);

load results52.mat;

thresholdVector = 0:0.005:0.2;
window_median_length_vector = 0;

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

save(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results' filesep 'pfa_data_median_64_lowT.mat'], 'averageFpr', 'stdFpr')

rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])
rmpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results']);  