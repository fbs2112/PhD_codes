clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '.' filesep 'Misc'])
if isunix
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep 'pfa_results' filesep]);
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep]);
    
else
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'statistical_data' filesep ...
        filesep 'pfa_results'])
    
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'statistical_data' filesep]);
end  

load results_det_2.mat;
load resultspfa1.mat;

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
dataPath = ['..' filesep '..' filesep '.' filesep 'figs' filesep 'group_presentation' filesep];

monteCarloLoops = 100;

thresholdVector = 0:0.05:0.5;

periodVector = 0;
bandwidthVector = 0;
JNRVector = 0;

for JNRIndex = 1:length(JNRVector)
    for bandwidthIndex = 1:length(bandwidthVector)
        for periodIndex = 1:length(periodVector)
            for thresholdIndex = 1:length(thresholdVector)
                for loopIndex = 1:monteCarloLoops
                    x = squeeze(detection_res(loopIndex, bandwidthIndex, periodIndex, JNRIndex, :, thresholdIndex));
                    tp(loopIndex, bandwidthIndex, periodIndex, JNRIndex, thresholdIndex) = median(x);
                    fn(loopIndex, bandwidthIndex, periodIndex, JNRIndex, thresholdIndex) = ...
                        1 - tp(loopIndex, bandwidthIndex, periodIndex, JNRIndex, thresholdIndex);
                end
            end
        end
    end
end

tpr = squeeze(tp./(tp+fn));

averageTpr = squeeze(mean(tpr, 1));
stdTpr = squeeze(std(tpr, [], 1));

cfun = @(tpr, fpr) sqrt(fpr.^2 + (1-tpr).^2);

for i = 1:length(JNRVector)
    c = cfun(averageTpr(i,:), averageFpr);
    [cMin(i), idx] = min(c(:));
    tprMin(i) = averageTpr(i,idx);
    fprMin(i) = averageFpr(idx);
end

figure;
plot(JNRVector, cMin)
ylim([0 1/sqrt(2)])
grid on;
ylabel('C$_{\mathrm{min}}$');
xlabel('JNR [dB]');

rmpath(['..' filesep '.' filesep 'Misc'])
if isunix
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep 'pfa_results' filesep]);
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep]);
    
else
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'statistical_data' filesep ...
        filesep 'pfa_results'])
    
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'statistical_data' filesep]);
end  