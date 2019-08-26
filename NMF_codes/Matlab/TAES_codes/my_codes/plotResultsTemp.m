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
addpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results']);  

load pfa_data_median_full_64_lowT_2;

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
dataPath = ['..' filesep '..' filesep '.' filesep 'figs' filesep 'group_presentation' filesep];

load results50
% end.mat;

monteCarloLoops = 100;

thresholdVector = 0:0.005:0.2;
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
rmpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results']); 


%-------------------------------Pai's results------------------------------

addpath(['..' filesep '..' filesep '.' filesep 'Misc'])
addpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'pai_results']);  

load pfa_data.mat;
load results11.mat;

averageFprPai = averageFpr;

fs = 32.768e6;

numberOfRawSamples = 4096;
silenceSamples = round(20e-6*fs);
totalSamples = numberOfRawSamples + silenceSamples*2;
WinLBlock = 19;
MBlock = fix(totalSamples/WinLBlock);

JNRVector = -25:0;

PfaVector = logspace(-5, 0, 17);
frequencySliceCW = 250;
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
    cMinPai(i) = min(c(:));
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
rmpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'pai_results']); 