clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '..' filesep '.' filesep 'Misc'])
addpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results']);  

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);

load results41.mat;

thresholdVector = 0.1:0.05:0.9;
window_median_length_vector = 51:50:401;
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
% for i = 1:size(averageFpr, 1)
plot(thresholdVector, averageFpr)
%     hold on;
% end

ylabel('Probability of false alarm');
xlabel('$\bar{\gamma}$');
xlim([min(thresholdVector) max(thresholdVector)]);

% legend('$L_{\mathrm{med}} = 51$', '$L_{\mathrm{med}} = 101$', '$L_{\mathrm{med}} = 151$','$L_{\mathrm{med}} = 201$',...
%     '$L_{\mathrm{med}} = 251$', '$L_{\mathrm{med}} = 301$', '$L_{\mathrm{med}} = 351$', '$L_{\mathrm{med}} = 401$');
grid on;

figure;
% for i = 1:size(averageFpr, 1)
    loglog(thresholdVector, averageFpr)
%     hold on;
% end

ylabel('Probability of false alarm');
xlabel('$\bar{\gamma}$');
xlim([min(thresholdVector) max(thresholdVector)]);

% legend('$L_{\mathrm{med}} = 51$', '$L_{\mathrm{med}} = 101$', '$L_{\mathrm{med}} = 151$','$L_{\mathrm{med}} = 201$',...
%     '$L_{\mathrm{med}} = 251$', '$L_{\mathrm{med}} = 301$', '$L_{\mathrm{med}} = 351$', '$L_{\mathrm{med}} = 401$');
grid on;

save(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results' filesep 'pfa_data_median_full_128_256.mat'], 'averageFpr', 'stdFpr')

rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])
rmpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results']);  