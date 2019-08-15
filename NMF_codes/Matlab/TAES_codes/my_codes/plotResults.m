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

% dataPath = ['..' filesep '..' filesep '.' filesep 'figs' filesep '08-12' filesep];

load results01.mat;

monteCarloLoops = 100;

thresholdVector = 0.1:0.05:0.9;
window_median_length_vector = 51:50:401;

for k = 1:length(thresholdVector)
    for j = 1:length(window_median_length_vector)
       
        x = squeeze(detection_res(:,k,j,:));
        fp(j,k,:) = sum(x, 2);
        tn(j,k,:) = size(x, 2) - fp(j,k,:) ;
    end
end

fpr = fp./(fp+tn);

averageFpr = squeeze(mean(fpr, 3));
stdFpr = squeeze(std(fpr, [], 3));


figure;
for i = 1:size(averageFpr, 1)
    plot(thresholdVector, averageFpr(i,:))
    hold on;
end

ylabel('Probability of false alarm');
xlabel('$\bar{\gamma}$');

legend('$L_{\mathrm{med}} = 51$', '$L_{\mathrm{med}} = 101$', '$L_{\mathrm{med}} = 151$','$L_{\mathrm{med}} = 201$',...
    '$L_{\mathrm{med}} = 251$', '$L_{\mathrm{med}} = 301$', '$L_{\mathrm{med}} = 351$', '$L_{\mathrm{med}} = 401$');
grid on;

figure;
for i = 1:size(averageFpr, 1)
    loglog(thresholdVector, averageFpr(i,:))
    hold on;
end

ylabel('Probability of false alarm');
xlabel('$\bar{\gamma}$');

legend('$L_{\mathrm{med}} = 51$', '$L_{\mathrm{med}} = 101$', '$L_{\mathrm{med}} = 151$','$L_{\mathrm{med}} = 201$',...
    '$L_{\mathrm{med}} = 251$', '$L_{\mathrm{med}} = 301$', '$L_{\mathrm{med}} = 351$', '$L_{\mathrm{med}} = 401$');
grid on;

rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])
rmpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results']);  