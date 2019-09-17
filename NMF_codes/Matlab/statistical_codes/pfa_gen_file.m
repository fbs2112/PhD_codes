clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '.' filesep 'Misc'])
if isunix
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
            'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
            filesep 'pfa_results' filesep]);
else
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'statistical_data' filesep ...
        filesep 'pfa_results' filesep])
end

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);

matNumber = 1;
load(['results' num2str(matNumber) '.mat']);

thresholdVector = 0:0.05:0.5;

for k = 1:length(thresholdVector)
    for i = 1:size(detection_res, 1)
        x = squeeze(detection_res(i,:,k));
        fp(i,k) = median(x);
        tn(i,k) = 1 - fp(i,k);
    end
end

fpr = fp./(fp+tn);
averageFpr = squeeze(mean(fpr, 1));
stdFpr = squeeze(std(fpr, [], 1));

for i = 1:length(thresholdVector)
    
    plot(thresholdVector, averageFpr)

    ylabel('Probability of false alarm');
    xlabel('$\bar{\gamma}$');
    xlim([min(thresholdVector) max(thresholdVector)]);
    grid on;
end

figure;
semilogy(thresholdVector, averageFpr)

ylabel('Probability of false alarm');
xlabel('$\bar{\gamma}$');
xlim([min(thresholdVector) max(thresholdVector)]);
grid on;

if isunix
    save(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
                'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
                filesep 'pfa_results' filesep 'resultspfa' num2str(matNumber) '.mat'], 'averageFpr', 'stdFpr')
else
    save(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'statistical_data' filesep ...
        filesep 'pfa_results' filesep 'resultspfa' num2str(matNumber) '.mat'], 'averageFpr', 'stdFpr')
end

rmpath(['..' filesep '.' filesep 'Misc'])
if isunix
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
            'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
            filesep 'pfa_results' filesep]);
else
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'statistical_data' filesep ...
        filesep 'pfa_results' filesep])
end