%%
%Plot results for CW considering the new approach of using a large median
%window

clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '..' filesep '.' filesep 'Misc'])
addpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results']);  

load pfa_data_median_full_512;

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
dataPath = ['..' filesep '..' filesep '.' filesep 'figs' filesep 'group_presentation' filesep];

load results28.mat;

monteCarloLoops = 100;

thresholdVector = 0.1:0.05:0.9;
window_median_length_vector = 0;
periodVector = 0;
bandwidthVector = 0;
JNRVector = -25:0;

for JNRIndex = 1:length(JNRVector)
    for bandwidthIndex = 1:length(bandwidthVector)
        for periodIndex = 1:length(periodVector)
            for thresholdIndex = 1:length(thresholdVector)
                x = squeeze(detection_res(:, bandwidthIndex, periodIndex, JNRIndex, thresholdIndex));
                tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :) = x;
                fn(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :) = ...
                    size(x, 2) - tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :);
            end
        end
    end
end

tpr = squeeze(tp./(tp+fn));

averageTpr = mean(tpr, 3);
stdTpr = std(tpr, [], 3);

JNRIdx = 1;
figure;
plot(thresholdVector, averageTpr(JNRIdx,:))

ylabel('Probability of detection');
xlabel('$\bar{\gamma}$');
xlim([min(thresholdVector) max(thresholdVector)]);
grid on;

figure;
loglog(thresholdVector, averageTpr(JNRIdx,:))

ylabel('Probability of detection');
xlabel('$\bar{\gamma}$');
xlim([min(thresholdVector) max(thresholdVector)]);
grid on;
    
for i = 1:size(averageTpr, 1)
    figure;
    plot(averageFpr, averageTpr(i,:));
    hold on 
    plot(linspace(0, 1, numel(averageFpr)), linspace(0, 1, numel(averageFpr)));
    title(['JNR: ' num2str(JNRVector(i))]);
    ylabel('Probability of detection');
    xlabel('Probability of false alarm');
    grid on;
    formatFig(gcf, [dataPath  'roc_CW_dot_' num2str(JNRVector(i))], 'en', figProp);
end

cfun = @(tpr, fpr) sqrt(fpr.^2 + (1-tpr).^2);

for i = 1:length(JNRVector)
    c = cfun(averageTpr(i,:), averageFpr);
    cMin(i) = min(c(:));
end


rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])
rmpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results']); 

%-------------------------------Pai's results------------------------------

addpath(['..' filesep '..' filesep '.' filesep 'Misc'])
addpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'pai_results']);  

load pfa_data.mat;
load results09.mat;

averageFprPai = averageFpr;

fs = 32.768e6;

numberOfRawSamples = 4096;
silenceSamples = round(20e-6*fs);
totalSamples = numberOfRawSamples + silenceSamples*2;
WinLBlock = 19;
MBlock = fix(totalSamples/WinLBlock);

JNRVector = -25:0;

PfaVector = logspace(-5, 0, 17);
frequencySliceCW = 28;
aux = 1:MBlock;
aux(frequencySliceCW) = 0;
aux = aux(aux~=0);

for JNRIndex = 1:length(JNRVector)
    
    for i = 1:length(PfaVector)
        x = squeeze(detection_res(JNRIndex,:,:,i));
        tpRegion = squeeze(x(:,frequencySliceCW,:));
        tpPai(JNRIndex,i,:) = tpRegion;
        fnPai(JNRIndex,i,:) = 1 - tpPai(JNRIndex,i,:);
    end

end

tprPai = tpPai./(tpPai+fnPai);
averageTprPai = mean(tprPai, 3);
stdTprPai = std(tprPai, [], 3);

for i = 1:size(averageTprPai, 1)
    figure;
    plot(averageFprPai(2,:), averageTprPai(i,:));
    hold on 
    plot(linspace(0, 1, numel(averageFpr)), linspace(0, 1, numel(averageFpr)));
    title(['JNR: ' num2str(JNRVector(i))]);
    ylabel('Probability of detection');
    xlabel('Probability of false alarm');
    grid on;
    formatFig(gcf, [dataPath  'roc_CW_pai_' num2str(JNRVector(i))], 'en', figProp);
end

for i = 1:length(JNRVector)
    c = cfun(averageTprPai(i,:), averageFprPai(2,:));
    cMinPai(i) = min(c(:));
end

figure;
plot(JNRVector, cMin)
hold on;
plot(JNRVector, cMinPai)

grid on;
ylabel('C$_{\mathrm{min}}$');
xlabel('JNR [dB]');
legend('Dot', 'Pai');

formatFig(gcf, [dataPath  'cmin_CW_'], 'en', figProp);

rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])
rmpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'pai_results']); 
%%
%Plot results for chirp considering the new approach of using a large median
%window

clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '..' filesep '.' filesep 'Misc'])
if isunix
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep 'pfa_results' filesep]);
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep]);
    
else
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep 'pfa_results' filesep])
    
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep]);
end  

load results_det_25.mat;
load resultspfa18.mat;

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
dataPath = ['..' filesep '..' filesep '.' filesep 'figs' filesep 'group_presentation' filesep];

monteCarloLoops = 100;

thresholdVector = 0:0.005:0.2;
% thresholdVector = 0.3:0.005:0.5;
% thresholdVector = -0.3:0.005:0;

periodVector = 0;
bandwidthVector = 0;
JNRVector = -25:0;

for JNRIndex = 1:length(JNRVector)
    for bandwidthIndex = 1:length(bandwidthVector)
        for periodIndex = 1:length(periodVector)
            for thresholdIndex = 1:length(thresholdVector)
                x = squeeze(detection_res(:, bandwidthIndex, periodIndex, JNRIndex, thresholdIndex));
                tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :) = x;
                fn(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :) = ...
                    size(x, 2) - tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :);
            end
        end
    end
end

tpr = squeeze(tp./(tp+fn));

averageTpr = mean(tpr, 3);
stdTpr = std(tpr, [], 3);

figure;
semilogy(thresholdVector, averageTpr(6:10,:))

ylabel('Probability of detection');
xlabel('$\bar{\gamma}$');
xlim([min(thresholdVector) max(thresholdVector)]);
grid on;


% for i = 1:size(averageTpr, 1)
%     figure;
%     plot(averageFpr, averageTpr(i,:));
%     hold on 
%     plot(linspace(0, 1, numel(averageFpr)), linspace(0, 1, numel(averageFpr)));
%     title(['JNR: ' num2str(JNRVector(i))]);
%     ylabel('Probability of detection');
%     xlabel('Probability of false alarm');
%     grid on;
% %     formatFig(gcf, [dataPath  'roc_chirp_dot_' num2str(JNRVector(i))], 'en', figProp)
% end

cfun = @(tpr, fpr) sqrt(fpr.^2 + (1-tpr).^2);

for i = 1:length(JNRVector)
    c = cfun(averageTpr(i,:), averageFpr);
    [cMin(i), idx] = min(c(:));
    tprMin(i) = averageTpr(i,idx);
    fprMin(i) = averageFpr(idx);
end

rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])
if isunix
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep 'pfa_results' filesep]);
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep]);
    
else
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep 'pfa_results' filesep])
    
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep]);
end  

%-------------------------------Pai's results------------------------------

addpath(['..' filesep '..' filesep '.' filesep 'Misc'])
if isunix
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep 'pfa_results' filesep]);
    
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep]);
    
else
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep 'pfa_results' filesep])
    
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep]);
end  

load resultspfa3.mat;
load results_det_5.mat;

averageFprPai = averageFpr;

JNRVector = -25:0;

PfaVector = logspace(-8, -2, 41);
% PfaVector = logspace(-5, -0, 17);

for JNRIndex = 1:length(JNRVector)
    
    for i = 1:length(PfaVector)
        tpRegion = squeeze(detection_res(JNRIndex,:,:,i));
        tpPai(JNRIndex,i,:) = any(tpRegion, 2);
        fnPai(JNRIndex,i,:) = 1 - tpPai(JNRIndex,i,:);
    end

end

tprPai = tpPai./(tpPai+fnPai);
averageTprPai = mean(tprPai, 3);
stdTprPai = std(tprPai, [], 3);

% for i = 1:size(averageTprPai, 1)
%     figure;
%     plot(averageFprPai(1,:), averageTprPai(i,:));
%     hold on 
%     plot(linspace(0, 1, numel(averageFpr)), linspace(0, 1, numel(averageFpr)));
%     title(['JNR: ' num2str(JNRVector(i))]);
%     ylabel('Probability of detection');
%     xlabel('Probability of false alarm');
%     grid on;
% %     formatFig(gcf, [dataPath  'roc_chirp_pai_' num2str(JNRVector(i))], 'en', figProp);
% end

for i = 1:length(JNRVector)
    c = cfun(averageTprPai(i,:), averageFprPai(1,:));
    [cMinPai(i), idx(i)] = min(c(:));
    tprMinPai(i) = averageTprPai(i,idx(i));
    fprMinPai(i) = averageFprPai(1,idx(i));
end


figure;
plot(JNRVector, cMin)
hold on;
plot(JNRVector, cMinPai)
ylim([0 1/sqrt(2)])
grid on;
ylabel('C$_{\mathrm{min}}$');
xlabel('JNR [dB]');
legend('Dot', 'Pai');

% formatFig(gcf, [dataPath  'cmin_chirp_'], 'en', figProp);


rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])
if isunix
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep 'pfa_results' filesep]);
    
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep]);
    
else
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep 'pfa_results' filesep])
    
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep]);
end  
%%
%Plot results for CW considering the new approach of using a large median
%window

clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '..' filesep '.' filesep 'Misc'])
if isunix
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep 'pfa_results' filesep]);
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep]);
    
else
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep 'pfa_results' filesep])
    
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep]);
end  

load results_det_24.mat;
load resultspfa18.mat;

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
dataPath = ['..' filesep '..' filesep '.' filesep 'figs' filesep 'group_presentation' filesep];

monteCarloLoops = 100;

thresholdVector = 0:0.005:0.2;
% thresholdVector = 0.3:0.005:0.5;
% thresholdVector = -0.3:0.005:0;

window_median_length_vector = 0;
periodVector = 0;
bandwidthVector = 0;
JNRVector = -25:0;

for JNRIndex = 1:length(JNRVector)
    for bandwidthIndex = 1:length(bandwidthVector)
        for periodIndex = 1:length(periodVector)
            for thresholdIndex = 1:length(thresholdVector)
                x = squeeze(detection_res(:, bandwidthIndex, periodIndex, JNRIndex, thresholdIndex));
                tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :) = x;
                fn(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :) = ...
                    size(x, 2) - tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :);
            end
        end
    end
end

tpr = squeeze(tp./(tp+fn));

averageTpr = mean(tpr, 3);
stdTpr = std(tpr, [], 3);

figure;
semilogy(thresholdVector, averageTpr(6:10,:))

ylabel('Probability of detection');
xlabel('$\bar{\gamma}$');
xlim([min(thresholdVector) max(thresholdVector)]);
grid on;


% for i = 1:size(averageTpr, 1)
%     figure;
%     plot(averageFpr, averageTpr(i,:));
%     hold on 
%     plot(linspace(0, 1, numel(averageFpr)), linspace(0, 1, numel(averageFpr)));
%     title(['JNR: ' num2str(JNRVector(i))]);
%     ylabel('Probability of detection');
%     xlabel('Probability of false alarm');
%     grid on;
% %     formatFig(gcf, [dataPath  'roc_chirp_dot_' num2str(JNRVector(i))], 'en', figProp)
% end

cfun = @(tpr, fpr) sqrt(fpr.^2 + (1-tpr).^2);

for i = 1:length(JNRVector)
    c = cfun(averageTpr(i,:), averageFpr);
    cMin(i) = min(c(:));
end

rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])
if isunix
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep 'pfa_results' filesep]);
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep]);
    
else
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep 'pfa_results' filesep])
    
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep]);
end  


%-------------------------------Pai's results------------------------------

addpath(['..' filesep '..' filesep '.' filesep 'Misc'])
if isunix
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep 'pfa_results' filesep]);
    
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep]);
    
else
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep 'pfa_results' filesep])
    
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep]);
end  

load resultspfa3.mat;
load results_det_7.mat;

averageFprPai = averageFpr;

JNRVector = -25:0;

PfaVector = logspace(-12, -2, 41);

for JNRIndex = 1:length(JNRVector)
    
    for i = 1:length(PfaVector)
        tpRegion = squeeze(detection_res(JNRIndex,:,:,i));
        tpPai(JNRIndex,i,:) = any(tpRegion, 2);
        fnPai(JNRIndex,i,:) = 1 - tpPai(JNRIndex,i,:);
    end

end

tprPai = tpPai./(tpPai+fnPai);
averageTprPai = mean(tprPai, 3);
stdTprPai = std(tprPai, [], 3);

% for i = 1:size(averageTprPai, 1)
%     figure;
%     plot(averageFprPai(1,:), averageTprPai(i,:));
%     hold on 
%     plot(linspace(0, 1, numel(averageFpr)), linspace(0, 1, numel(averageFpr)));
%     title(['JNR: ' num2str(JNRVector(i))]);
%     ylabel('Probability of detection');
%     xlabel('Probability of false alarm');
%     grid on;
% %     formatFig(gcf, [dataPath  'roc_chirp_pai_' num2str(JNRVector(i))], 'en', figProp);
% end

for i = 1:length(JNRVector)
    c = cfun(averageTprPai(i,:), averageFpr(1,:));
    [cMinPai(i), idx(i)] = min(c(:));
end


figure;
plot(JNRVector, cMin)
hold on;
plot(JNRVector, cMinPai)
ylim([0 1/sqrt(2)])
grid on;
ylabel('C$_{\mathrm{min}}$');
xlabel('JNR [dB]');
legend('Dot', 'Pai');

% formatFig(gcf, [dataPath  'cmin_chirp_'], 'en', figProp);


rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])
if isunix
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep 'pfa_results' filesep]);
    
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep]);
    
else
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep 'pfa_results' filesep])
    
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep]);
end  

%%
%Plot results for chirp for different bandwidths and periods

clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '..' filesep '.' filesep 'Misc'])
if isunix
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep 'pfa_results' filesep]);
    
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep]);
    
else
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep 'pfa_results' filesep])
    
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep]);
end  

load results_det_23.mat;
load resultspfa14.mat;

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
dataPath = ['..' filesep '..' filesep '.' filesep 'figs' filesep 'group_presentation' filesep];

monteCarloLoops = 100;

thresholdVector = 0:0.005:0.2;

bandwidthVector = (2e6:3e6:14e6)/1e6;
periodVector = (8.62e-6:1.48e-6:18.97e-6)*1e6;

JNRVector = -25:0;

for JNRIndex = 1:length(JNRVector)
    for bandwidthIndex = 1:length(bandwidthVector)
        for periodIndex = 1:length(periodVector)
            for thresholdIndex = 1:length(thresholdVector)
                x = squeeze(detection_res(:, bandwidthIndex, periodIndex, JNRIndex, thresholdIndex));
                tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :) = x;
                fn(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :) = ...
                    size(x, 2) - tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :);
            end
        end
    end
end

tpr = squeeze(tp./(tp+fn));

averageTpr = mean(tpr, 5);
stdTpr = std(tpr, [], 5);

cfun = @(tpr, fpr) sqrt(fpr.^2 + (1-tpr).^2);

for k = 1:length(bandwidthVector)
    for j = 1:length(periodVector)
        for i = 1:length(JNRVector)
            c = cfun(squeeze(averageTpr(k,j,i,:)).', averageFpr);
            [cMin(k,j,i), idx] = min(c(:));
            tprMin(k,j,i) = averageTpr(k,j,i,idx);
            fprMin(k,j,i) = averageFpr(idx);
        end
    end
end

for k = 1:length(bandwidthVector)
    figure;
    
    for j = 1:length(periodVector)
        plot(JNRVector, squeeze(cMin(k,j,:)));
        title([num2str(bandwidthVector(k)) ' MHz'])
        hold on;
    end
    ylabel('C$_{\mathrm{min}}$');
    xlabel('JNR [dB]');
    legend(['T = ' num2str(periodVector(1)) ' $\mu$s'], ['T = ' num2str(periodVector(2)) ' $\mu$s'], ['T = ' num2str(periodVector(3)) ' $\mu$s'],...
        ['T = ' num2str(periodVector(4)) ' $\mu$s'], ['T = ' num2str(periodVector(5)) ' $\mu$s'], ['T = ' num2str(periodVector(6)) ' $\mu$s'],...
        ['T = ' num2str(periodVector(7)) ' $\mu$s']);
    grid on;
end

rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])
if isunix
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep 'pfa_results' filesep]);
    
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep]);
    
else
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep 'pfa_results' filesep])
    
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep]);
end  

%%
%Plot results for chirp for different bandwidths and periods for Pai's
%technique

clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '..' filesep '.' filesep 'Misc'])
if isunix
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep 'pfa_results' filesep]);
    
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep]);
    
else
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep 'pfa_results' filesep])
    
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep]);
end  

load results_det_10.mat;
load resultspfa3.mat;

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
dataPath = ['..' filesep '..' filesep '.' filesep 'figs' filesep 'group_presentation' filesep];

monteCarloLoops = 100;

PfaVector = logspace(-12, -2, 41);

bandwidthVector = (2e6:3e6:14e6)/1e6;
periodVector = (8.62e-6:1.48e-6:18.97e-6)*1e6;

JNRVector = -25:0;

for JNRIndex = 1:length(JNRVector)
    for bandwidthIndex = 1:length(bandwidthVector)
        for periodIndex = 1:length(periodVector)
            for thresholdIndex = 1:length(PfaVector)
                x = squeeze(detection_res(bandwidthIndex, periodIndex, JNRIndex, :, :, thresholdIndex));
                tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :) = any(x, 2);
                fn(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :) = ...
                    1 - any(x, 2);
            end
        end
    end
end

tpr = squeeze(tp./(tp+fn));

averageTpr = mean(tpr, 5);
stdTpr = std(tpr, [], 5);

cfun = @(tpr, fpr) sqrt(fpr.^2 + (1-tpr).^2);

for k = 1:length(bandwidthVector)
    for j = 1:length(periodVector)
        for i = 1:length(JNRVector)
            c = cfun(squeeze(averageTpr(k,j,i,:)).', averageFpr(:,2));
            [cMin(k,j,i), idx] = min(c(:));
        end
    end
end

for k = 1:length(bandwidthVector)
    figure;
    
    for j = 1:length(periodVector)
        plot(JNRVector, squeeze(cMin(k,j,:)));
        title([num2str(bandwidthVector(k)) ' MHz'])
        hold on;
    end
    ylabel('C$_{\mathrm{min}}$');
    xlabel('JNR [dB]');
    legend(['T = ' num2str(periodVector(1)) ' $\mu$s'], ['T = ' num2str(periodVector(2)) ' $\mu$s'], ['T = ' num2str(periodVector(3)) ' $\mu$s'],...
        ['T = ' num2str(periodVector(4)) ' $\mu$s'], ['T = ' num2str(periodVector(5)) ' $\mu$s'], ['T = ' num2str(periodVector(6)) ' $\mu$s'],...
        ['T = ' num2str(periodVector(7)) ' $\mu$s']);
    grid on;
end

rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])

if isunix
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep 'pfa_results' filesep]);
    
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep]);
    
else
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep 'pfa_results' filesep])
    
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep]);
end  