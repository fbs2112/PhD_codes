clear;
clc;
close all;

addpath(['.' filesep 'data']);

load resultsPai01.mat

PfaVector = logspace(-5, 0, 25);

for i = 1:length(PfaVector)
   x = squeeze(GoFBlockDeteflag(i,:,:));
   tpRegion = x(26,:);
   fpRegion = x(50:end,:);
   for j = 1:1
        tp(j,i) = tpRegion(j);
        fn(j,i) = 1 - tp(j,i);
        fp(j,i) = sum(fpRegion(:,j));
        tn(j,i) = size(fpRegion(:,j), 1) - fp(j,i);
   end
end

tpr = tp./(tp+fn);
fpr = fp./(fp+tn);

tprAverage = mean(tpr, 1);
fprAverage = mean(fpr, 1);
tprStd = std(tpr, 1);
fprStd = std(fpr, 1);


figure;
plot(fprAverage, tprAverage);
hold on;
xAux = linspace(0, 1, numel(fprAverage));
yAux = linspace(0, 1, numel(fprAverage));
plot(xAux, yAux, '--');
box 'off'
grid on
xlabel('Probability of false alarm');
ylabel('Probability of detection');

figure;
loglog(fprAverage, tprAverage);
xlabel('Probability of false alarm');
ylabel('Probability of detection');

figure;
loglog(PfaVector, tprAverage);
xlabel('Probability of false alarm');
ylabel('Probability of detection');

rmpath(['.' filesep 'data']);
%%

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

clear;
clc;
close all;

addpath(['..' filesep '..' filesep '.' filesep 'Misc'])
addpath(['..' filesep '..' filesep '.' filesep 'data' filesep '08-12']);  

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);

dataPath = ['..' filesep '..' filesep '.' filesep 'figs' filesep '08-12' filesep];

load results04.mat;

onset = 528;
offset = 4623;
window_length = 98;

thresholdVector = 0.1:0.05:0.9;
window_median_length_vector = 51:50:401;
monteCarloLoops = 100;
JNRVector = -17;

bandwidthVector = 10.72e6;
periodVector = 8.72e-6;

averageFP = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
averageTN = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
averageTP = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
averageFN = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));

averageTPR = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
averageFPR = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
stdTPR = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
stdFPR = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));

averageACC = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
stdACC = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
averageF1 = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
stdF1 = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
averagePrec = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
stdPrec = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
averageCmin = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
stdCmin = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));


cminfun = @(tpr, fpr) sqrt(fpr.^2 + (1-tpr).^2);

for JNRIndex = 1:length(JNRVector)
    for bandwidthVectorIndex = 1:length(bandwidthVector)
        for periodIndex = 1:length(periodVector)
            for thresholdIndex = 1:length(thresholdVector)
                for window_median_length_index = 1:length(window_median_length_vector)
                    
                    fp = zeros(monteCarloLoops, 1);
                    tp = zeros(monteCarloLoops, 1);
                    fn = zeros(monteCarloLoops, 1);
                    tn = zeros(monteCarloLoops, 1);
                    
                    for loopIndex = 1:monteCarloLoops
                        fpRegion1 = squeeze(detection_res(loopIndex, bandwidthVectorIndex, periodIndex, JNRIndex, thresholdIndex, window_median_length_index, 1:onset));
                        tpRegion = squeeze(detection_res(loopIndex, bandwidthVectorIndex, periodIndex, JNRIndex, thresholdIndex, window_median_length_index, onset  + round(window_length*2):offset - round(window_length*2)));
                        fpRegion2 = squeeze(detection_res(loopIndex, bandwidthVectorIndex, periodIndex, JNRIndex, thresholdIndex, window_median_length_index, offset + round(window_length*2):end));
                        fpRegion = [fpRegion1;fpRegion2];
                        
                        fp(loopIndex) = sum(fpRegion);
                        tn(loopIndex) = length(fpRegion) - fp(loopIndex);
                        
                        tp(loopIndex) = sum(tpRegion);
                        fn(loopIndex) = length(tpRegion) - tp(loopIndex);
                    end
                    
                    acc = (tp + tn)./(tp + fn + tn + fp);
                    f1Score = 2*tp./(2*tp + fp + fn);
                    precision = tp./(tp+fp);
                    precision(isnan(precision)) = 0;
                    tpr = tp./(tp+fn);
                    fpr = fp./(fp+tn);
                    
                    
                    cmin = cminfun(tpr, fpr);
                    
                    fpStd = std(fp);
                    tnStd = std(tn);
                    tpStd = std(tp);
                    fnStd = std(fn);
                    
                    averageFP(JNRIndex, thresholdIndex, window_median_length_index) = mean(fp);
                    averageTN(JNRIndex, thresholdIndex, window_median_length_index) = mean(tn);
                    averageTP(JNRIndex, thresholdIndex, window_median_length_index) = mean(tp);
                    averageFN(JNRIndex, thresholdIndex, window_median_length_index) = mean(fn);
                    
                    averageTPR(JNRIndex, thresholdIndex, window_median_length_index) = mean(tpr);
                    averageFPR(JNRIndex, thresholdIndex, window_median_length_index) = mean(fpr);
                    stdTPR(JNRIndex, thresholdIndex, window_median_length_index) = std(tpr);
                    stdFPR(JNRIndex, thresholdIndex, window_median_length_index) = std(fpr);
                    
                    averageACC(JNRIndex, thresholdIndex, window_median_length_index) = mean(acc);
                    stdACC(JNRIndex, thresholdIndex, window_median_length_index) = std(acc);
                    averageF1(JNRIndex, thresholdIndex, window_median_length_index) = mean(f1Score);
                    stdF1(JNRIndex, thresholdIndex, window_median_length_index) = std(f1Score);
                    averagePrec(JNRIndex, thresholdIndex, window_median_length_index) = mean(precision);
                    stdPrec(JNRIndex, thresholdIndex, window_median_length_index) = std(precision);
                    
                    averageCmin(JNRIndex, thresholdIndex, window_median_length_index) = mean(cmin);
                    stdCmin(JNRIndex, thresholdIndex, window_median_length_index) = std(cmin);
                end
            end
        end
    end
end

figure;
for i = 1:length(JNRVector)
    averageACCAux = squeeze(averageACC(i,:,:));

    [maxval, maxidx] = max(averageACCAux(:));
    [j, k] = ind2sub(size(averageACCAux), maxidx);
    disp('---------------------Accuracy-----------------------------------')
    disp('------------------------------------')
    disp(['-------------JNR: ' num2str(JNRVector(i)) ' dB-----------------------'])
    disp(['Accuracy: ' num2str(maxval) '+-' num2str(stdACC(i,j,k))])
    disp('------------------------------------')
    disp(['Pd: ' num2str(averageTPR(i,j,k)) '+-' num2str(stdTPR(i,j,k))])
    disp('------------------------------------')
    disp(['Pfa: ' num2str(averageFPR(i,j,k)) '+-' num2str(stdFPR(i,j,k))])
    disp('------------------------------------')
    disp(['Threshold: ' num2str(thresholdVector(j))])
    disp('------------------------------------')
    disp(['Window length: ' num2str(window_median_length_vector(k))])
    disp('------------------------------------')
    
    plot(averageFPR(i,:,k), averageTPR(i,:,k))
    hold on    
end

for i = 1:length(JNRVector)
    averageF1Aux = squeeze(averageF1(i,:,:));

    [maxval, maxidx] = max(averageF1Aux(:));
    [j, k] = ind2sub(size(averageF1Aux), maxidx);
    disp('---------------------F1-score-----------------------------------')
    disp('------------------------------------')
    disp(['-------------JNR: ' num2str(JNRVector(i)) ' dB-----------------------'])
    disp(['F1-score: ' num2str(maxval) '+-' num2str(stdF1(i,j,k))])
    disp('------------------------------------')
    disp(['Pd: ' num2str(averageTPR(i,j,k)) '+-' num2str(stdTPR(i,j,k))])
    disp('------------------------------------')
    disp(['Pfa: ' num2str(averageFPR(i,j,k)) '+-' num2str(stdFPR(i,j,k))])
    disp('------------------------------------')
    disp(['Threshold: ' num2str(thresholdVector(j))])
    disp('------------------------------------')
    disp(['Window length: ' num2str(window_median_length_vector(k))])
    disp('------------------------------------')
    
    plot(averageFPR(i,:,k), averageTPR(i,:,k))
end


for i = 1:length(JNRVector)
    averagePrecAux = squeeze(averagePrec(i,:,:));

    [maxval, maxidx] = max(averagePrecAux(:), [], 'includenan');
    [j, k] = ind2sub(size(averagePrecAux), maxidx);
    disp('---------------------Precision-----------------------------------')
    disp('------------------------------------')
    disp(['-------------JNR: ' num2str(JNRVector(i)) ' dB-----------------------'])
    disp(['Precision: ' num2str(maxval) '+-' num2str(stdPrec(i,j,k))])
    disp('------------------------------------')
    disp(['Pd: ' num2str(averageTPR(i,j,k)) '+-' num2str(stdTPR(i,j,k))])
    disp('------------------------------------')
    disp(['Pfa: ' num2str(averageFPR(i,j,k)) '+-' num2str(stdFPR(i,j,k))])
    disp('------------------------------------')
    disp(['Threshold: ' num2str(thresholdVector(j))])
    disp('------------------------------------')
    disp(['Window length: ' num2str(window_median_length_vector(k))])
    disp('------------------------------------')
    
    plot(averageFPR(i,:,k), averageTPR(i,:,k))
end

for i = 1:length(JNRVector)

    averageCminAux = squeeze(averageCmin(i,:,:));
    [minCmin, minIdx] = min(averageCminAux(:));
    [j, k] = ind2sub(size(averageCminAux), minIdx);

    disp('---------------------Cmin-----------------------------------')
    disp('------------------------------------')
    disp(['-------------JNR: ' num2str(JNRVector(i)) ' dB-----------------------'])
    disp(['Cmin: ' num2str(minCmin) '+-' num2str(stdCmin(i,j,k))])
    disp('------------------------------------')
    disp(['Pd: ' num2str(averageTPR(i,j,k)) '+-' num2str(stdTPR(i,j,k))])
    disp('------------------------------------')
    disp(['Pfa: ' num2str(averageFPR(i,j,k)) '+-' num2str(stdFPR(i,j,k))])
    disp('------------------------------------')
    disp(['Threshold: ' num2str(thresholdVector(j))])
    disp('------------------------------------')
    disp(['Window length: ' num2str(window_median_length_vector(k))])
    disp('------------------------------------')
end

clc;

addpath(['.' filesep 'data']);

% load resultsPai02_2.mat
load resultsPai01.mat


PfaVector = logspace(-5, 0, 25);
frequencySliceCW = 26;
% frequencySliceCW = 35;

aux = 1:size(GoFBlockDeteflag, 2);
aux(frequencySliceCW) = 0;
aux = aux(aux~=0);

for i = 1:length(PfaVector)
   x = squeeze(GoFBlockDeteflag(i,:,:));
   tpRegion = x(frequencySliceCW,:);
   fpRegion = x(aux,:);
   for j = 1:size(x, 2)
       tp(j,i) = tpRegion(j);
       fn(j,i) = 1 - tp(j,i);
       fp(j,i) = sum(fpRegion(:,j));
       tn(j,i) = size(fpRegion(:,j), 1) - fp(j,i);
   end
end

tpr = tp./(tp+fn);
fpr = fp./(fp+tn);

acc = (tp + tn)./(tp + fn + tn + fp);
f1Score = 2*tp./(2*tp + fp + fn);
precision = tp./(tp+fp);
precision(isnan(precision)) = 0;

tpAverage = mean(tp, 1);
tnAverage = mean(tn, 1);
fpAverage = mean(fp, 1);
fnAverage = mean(fn, 1);

tprAverage = mean(tpr, 1);
fprAverage = mean(fpr, 1);
tprStd = std(tpr, [], 1);
fprStd = std(fpr, [], 1);

accAverage = mean(acc, 1);
f1ScoreAverage = mean(f1Score, 1);
precAverage = mean(precision, 1);
accStd = std(acc, [], 1);
f1ScoreStd = std(f1Score, [], 1);
precStd = std(precision, [], 1);

% figure;
% errorbar(fprAverage,tprAverage, fprStd, tprStd, 'both');
plot(fprAverage, tprAverage);
% hold on;

% loglog(PfaVector, tprAverage);
% hold on;

xAux = linspace(0, 1, numel(fprAverage));
yAux = linspace(0, 1, numel(fprAverage));
plot(xAux, yAux, '--');

% xAux = logspace(-5, 0, numel(fprAverage));
% yAux = logspace(-5, 0, numel(fprAverage));
% loglog(xAux, yAux, '--');
box 'off'
% xlim([1e-5 1e0])
% ylim([1e-5 1e0])

grid on
xlabel('Probability of false alarm');
% xlabel('$\bar{\gamma}$');

ylabel('Probability of detection');
legend('Acc.', 'F$_1$-score', 'Prec.', 'Statistical', 'Random guess', 'Location', 'southeast');
% legend('Statistical', 'Random guess', 'Location', 'southeast');

[maxval, maxidx] = max(accAverage);
disp('---------------------Accuracy-----------------------------------')
disp('------------------------------------')
disp(['Accuracy: ' num2str(maxval) '+-' num2str(accStd(maxidx))])
disp('------------------------------------')
disp(['Pd: ' num2str(tprAverage(maxidx)) '+-' num2str(tprStd(maxidx))])
disp('------------------------------------')
disp(['Pfa: ' num2str(fprAverage(maxidx)) '+-' num2str(fprStd(maxidx))])
disp('------------------------------------')
disp(['Threshold: ' num2str(PfaVector(maxidx))])
disp('------------------------------------')

[maxval, maxidx] = max(precAverage);
disp('---------------------Precision-----------------------------------')
disp('------------------------------------')
disp(['Precision: ' num2str(maxval) '+-' num2str(precStd(maxidx))])
disp('------------------------------------')
disp(['Pd: ' num2str(tprAverage(maxidx)) '+-' num2str(tprStd(maxidx))])
disp('------------------------------------')
disp(['Pfa: ' num2str(fprAverage(maxidx)) '+-' num2str(fprStd(maxidx))])
disp('------------------------------------')
disp(['Threshold: ' num2str(PfaVector(maxidx))])
disp('------------------------------------')

[maxval, maxidx] = max(tprAverage);
disp('---------------------Recall-----------------------------------')
disp('------------------------------------')
disp(['Recall: ' num2str(maxval) '+-' num2str(tprStd(maxidx))])
disp('------------------------------------')
disp(['Pd: ' num2str(tprAverage(maxidx)) '+-' num2str(tprStd(maxidx))])
disp('------------------------------------')
disp(['Pfa: ' num2str(fprAverage(maxidx)) '+-' num2str(fprStd(maxidx))])
disp('------------------------------------')
disp(['Threshold: ' num2str(PfaVector(maxidx))])
disp('------------------------------------')

[maxval, maxidx] = max(f1ScoreAverage);
disp('---------------------F1-score-----------------------------------')
disp('------------------------------------')
disp(['F1-score: ' num2str(maxval) '+-' num2str(f1ScoreStd(maxidx))])
disp('------------------------------------')
disp(['Pd: ' num2str(tprAverage(maxidx)) '+-' num2str(tprStd(maxidx))])
disp('------------------------------------')
disp(['Pfa: ' num2str(fprAverage(maxidx)) '+-' num2str(fprStd(maxidx))])
disp('------------------------------------')
disp(['Threshold: ' num2str(PfaVector(maxidx))])
disp('------------------------------------')


cminfun = @(tpr, fpr) sqrt(fpr.^2 + (1-tpr).^2);

cmin = cminfun(tpr, fpr);

cminAverage = mean(cmin, 1);
cminStd = std(cmin, [], 1);

[minCmin, minIdx] = min(cminAverage);

disp('---------------------Cmin-----------------------------------')
disp('------------------------------------')
disp(['Cmin: ' num2str(minCmin) '+-' num2str(cminStd(minIdx))])
disp('------------------------------------')
disp(['Pd: ' num2str(tprAverage(minIdx)) '+-' num2str(tprStd(minIdx))])
disp('------------------------------------')
disp(['Pfa: ' num2str(fprAverage(minIdx)) '+-' num2str(fprStd(minIdx))])
disp('------------------------------------')
disp(['Threshold: ' num2str(PfaVector(minIdx))])
disp('------------------------------------')
formatFig(gcf, [dataPath  'roc_CW_2'], 'en', figProp);


rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])
rmpath(['..' filesep '..' filesep '.' filesep 'data' filesep '08-12']);  
rmpath(['.' filesep 'data']);

%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

clear;
clc;
close all;

addpath(['..' filesep '..' filesep '.' filesep 'Misc'])
addpath(['..' filesep '..' filesep '.' filesep 'data' filesep '08-12']);  

onset = 528;
offset = 4623;
window_length = 98;

thresholdVector = 0.1:0.05:0.9;
window_median_length_vector = 51:50:401;
monteCarloLoops = 100;
JNRVector = -20:0;

bandwidthVector = 10.72e6;
periodVector = 8.72e-6;

averageFP = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
averageTN = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
averageTP = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
averageFN = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));

averageTPR = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
averageFPR = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
stdTPR = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
stdFPR = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));

averageACC = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
stdACC = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
averageF1 = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
stdF1 = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
averagePrec = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
stdPrec = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
averageCmin = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
stdCmin = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));

cminfun = @(tpr, fpr) sqrt(fpr.^2 + (1-tpr).^2);

for JNRIndex = 1:length(JNRVector)
    load(['results07_' num2str(JNRIndex)]);
    for bandwidthVectorIndex = 1:length(bandwidthVector)
        for periodIndex = 1:length(periodVector)
            for thresholdIndex = 1:length(thresholdVector)
                for window_median_length_index = 1:length(window_median_length_vector)
                    
                    fp = zeros(monteCarloLoops, 1);
                    tp = zeros(monteCarloLoops, 1);
                    fn = zeros(monteCarloLoops, 1);
                    tn = zeros(monteCarloLoops, 1);
                    
                    for loopIndex = 1:monteCarloLoops
                        fpRegion1 = squeeze(detection_res(loopIndex, bandwidthVectorIndex, periodIndex, 1, thresholdIndex, window_median_length_index, 1:onset));
                        tpRegion = squeeze(detection_res(loopIndex, bandwidthVectorIndex, periodIndex, 1, thresholdIndex, window_median_length_index, onset  + round(window_length*2):offset - round(window_length*2)));
                        fpRegion2 = squeeze(detection_res(loopIndex, bandwidthVectorIndex, periodIndex, 1, thresholdIndex, window_median_length_index, offset + round(window_length*2):end));
                        fpRegion = [fpRegion1;fpRegion2];
                        
                        fp(loopIndex) = sum(fpRegion);
                        tn(loopIndex) = length(fpRegion) - fp(loopIndex);
                        
                        tp(loopIndex) = sum(tpRegion);
                        fn(loopIndex) = length(tpRegion) - tp(loopIndex);
                    end
                    
                    acc = (tp + tn)./(tp + fn + tn + fp);
                    f1Score = 2*tp./(2*tp + fp + fn);
                    precision = tp./(tp+fp);
                    precision(isnan(precision)) = 0;
                    tpr = tp./(tp+fn);
                    fpr = fp./(fp+tn);
                    
                    cmin = cminfun(tpr, fpr);
                    
                    fpStd = std(fp);
                    tnStd = std(tn);
                    tpStd = std(tp);
                    fnStd = std(fn);
                    
                    averageFP(JNRIndex, thresholdIndex, window_median_length_index) = mean(fp);
                    averageTN(JNRIndex, thresholdIndex, window_median_length_index) = mean(tn);
                    averageTP(JNRIndex, thresholdIndex, window_median_length_index) = mean(tp);
                    averageFN(JNRIndex, thresholdIndex, window_median_length_index) = mean(fn);
                    
                    averageTPR(JNRIndex, thresholdIndex, window_median_length_index) = mean(tpr);
                    averageFPR(JNRIndex, thresholdIndex, window_median_length_index) = mean(fpr);
                    stdTPR(JNRIndex, thresholdIndex, window_median_length_index) = std(tpr);
                    stdFPR(JNRIndex, thresholdIndex, window_median_length_index) = std(fpr);
                    
                    averageACC(JNRIndex, thresholdIndex, window_median_length_index) = mean(acc);
                    stdACC(JNRIndex, thresholdIndex, window_median_length_index) = std(acc);
                    averageF1(JNRIndex, thresholdIndex, window_median_length_index) = mean(f1Score);
                    stdF1(JNRIndex, thresholdIndex, window_median_length_index) = std(f1Score);
                    averagePrec(JNRIndex, thresholdIndex, window_median_length_index) = mean(precision);
                    stdPrec(JNRIndex, thresholdIndex, window_median_length_index) = std(precision);
                    
                    averageCmin(JNRIndex, thresholdIndex, window_median_length_index) = mean(cmin);
                    stdCmin(JNRIndex, thresholdIndex, window_median_length_index) = std(cmin);
                end
            end
        end
    end
end

for i = 1:length(JNRVector)
    averageACCAux = squeeze(averageACC(i,:,:));

    [maxval, maxidx] = max(averageACCAux(:));
    [j, k] = ind2sub(size(averageACCAux), maxidx);
    disp('---------------------Accuracy-----------------------------------')
    disp('------------------------------------')
    disp(['-------------JNR: ' num2str(JNRVector(i)) ' dB-----------------------'])
    disp(['Accuracy: ' num2str(maxval) '+-' num2str(stdACC(i,j,k))])
    disp('------------------------------------')
    disp(['Pd: ' num2str(averageTPR(i,j,k)) '+-' num2str(stdTPR(i,j,k))])
    disp('------------------------------------')
    disp(['Pfa: ' num2str(averageFPR(i,j,k)) '+-' num2str(stdFPR(i,j,k))])
    disp('------------------------------------')
    disp(['Threshold: ' num2str(thresholdVector(j))])
    disp('------------------------------------')
    disp(['Window length: ' num2str(window_median_length_vector(k))])
    disp('------------------------------------')
    
%     plot(averageFPR(i,:,k), averageTPR(i,:,k))
    hold on    
end

for i = 1:length(JNRVector)
    averageF1Aux = squeeze(averageF1(i,:,:));

    [maxval, maxidx] = max(averageF1Aux(:));
    [j, k] = ind2sub(size(averageF1Aux), maxidx);
    disp('---------------------F1-score-----------------------------------')
    disp('------------------------------------')
    disp(['-------------JNR: ' num2str(JNRVector(i)) ' dB-----------------------'])
    disp(['F1-score: ' num2str(maxval) '+-' num2str(stdF1(i,j,k))])
    disp('------------------------------------')
    disp(['Pd: ' num2str(averageTPR(i,j,k)) '+-' num2str(stdTPR(i,j,k))])
    disp('------------------------------------')
    disp(['Pfa: ' num2str(averageFPR(i,j,k)) '+-' num2str(stdFPR(i,j,k))])
    disp('------------------------------------')
    disp(['Threshold: ' num2str(thresholdVector(j))])
    disp('------------------------------------')
    disp(['Window length: ' num2str(window_median_length_vector(k))])
    disp('------------------------------------')
    
%     plot(averageFPR(i,:,k), averageTPR(i,:,k))
end


for i = 1:length(JNRVector)
    averagePrecAux = squeeze(averagePrec(i,:,:));

    [maxval, maxidx] = max(averagePrecAux(:), [], 'includenan');
    [j, k] = ind2sub(size(averagePrecAux), maxidx);
    disp('---------------------Precision-----------------------------------')
    disp('------------------------------------')
    disp(['-------------JNR: ' num2str(JNRVector(i)) ' dB-----------------------'])
    disp(['Precision: ' num2str(maxval) '+-' num2str(stdPrec(i,j,k))])
    disp('------------------------------------')
    disp(['Pd: ' num2str(averageTPR(i,j,k)) '+-' num2str(stdTPR(i,j,k))])
    disp('------------------------------------')
    disp(['Pfa: ' num2str(averageFPR(i,j,k)) '+-' num2str(stdFPR(i,j,k))])
    disp('------------------------------------')
    disp(['Threshold: ' num2str(thresholdVector(j))])
    disp('------------------------------------')
    disp(['Window length: ' num2str(window_median_length_vector(k))])
    disp('------------------------------------')
    
%     plot(averageFPR(i,:,k), averageTPR(i,:,k))
end

for i = 1:length(JNRVector)

    averageCminAux = squeeze(averageCmin(i,:,:));
    [minCminNMF(i), minIdx] = min(averageCminAux(:));
    [j, k] = ind2sub(size(averageCminAux), minIdx);

    disp('---------------------Cmin-----------------------------------')
    disp('------------------------------------')
    disp(['-------------JNR: ' num2str(JNRVector(i)) ' dB-----------------------'])
    disp(['Cmin: ' num2str(minCminNMF(i)) '+-' num2str(stdCmin(i,j,k))])
    disp('------------------------------------')
    disp(['Pd: ' num2str(averageTPR(i,j,k)) '+-' num2str(stdTPR(i,j,k))])
    disp('------------------------------------')
    disp(['Pfa: ' num2str(averageFPR(i,j,k)) '+-' num2str(stdFPR(i,j,k))])
    disp('------------------------------------')
    disp(['Threshold: ' num2str(thresholdVector(j))])
    disp('------------------------------------')
    disp(['Window length: ' num2str(window_median_length_vector(k))])
    disp('------------------------------------')
end

clc;

addpath(['.' filesep 'data']);

fs = 32.768e6;

numberOfRawSamples = 4096;
silenceSamples = round(20e-6*fs);
totalSamples = numberOfRawSamples + silenceSamples*2;
WinLBlock = 19;
MBlock = fix(totalSamples/WinLBlock);

JNRVector = -20:0;

PfaVector = logspace(-5, 0, 17);
frequencySliceCW = 35;
aux = 1:MBlock;
aux(frequencySliceCW) = 0;
aux = aux(aux~=0);

for JNRIndex = 1:length(JNRVector)
    
    load(['resultsPai03_' num2str(JNRIndex) '.mat']);
    for i = 1:length(PfaVector)
        x = squeeze(GoFBlockDeteflag(i,:,:));
        tpRegion = x(frequencySliceCW,:);
        fpRegion = x(aux,:);
        for j = 1:size(x, 2)
            tp(JNRIndex,j,i) = tpRegion(j);
            fn(JNRIndex,j,i) = 1 - tp(JNRIndex,j,i);
            fp(JNRIndex,j,i) = sum(fpRegion(:,j));
            tn(JNRIndex,j,i) = size(fpRegion(:,j), 1) - fp(JNRIndex,j,i);
        end
    end

end

tpr = tp./(tp+fn);
fpr = fp./(fp+tn);

tprAverage = squeeze(mean(tpr, 2));
fprAverage = squeeze(mean(fpr, 2));
tprStd = squeeze(std(tpr, [], 2));
fprStd = squeeze(std(fpr, [], 2));

cminfun = @(tpr, fpr) sqrt(fpr.^2 + (1-tpr).^2);
cmin = cminfun(tpr, fpr);
cminAverage = squeeze(mean(cmin, 2));
cminStd = squeeze(std(cmin, [], 2));    
    

for i = 1:length(JNRVector)
    [minCmin(i), minIdx] = min(cminAverage(i,:));

    disp('---------------------Cmin-----------------------------------')
    disp(['JNR: ' num2str(JNRVector(i)) ' dB'])
    disp('------------------------------------')
    disp(['Cmin: ' num2str(minCmin(i)) '+-' num2str(cminStd(i,minIdx))])
    disp('------------------------------------')
    disp(['Pd: ' num2str(tprAverage(i,minIdx)) '+-' num2str(tprStd(i,minIdx))])
    disp('------------------------------------')
    disp(['Pfa: ' num2str(fprAverage(i,minIdx)) '+-' num2str(fprStd(i,minIdx))])
    disp('------------------------------------')
    disp(['Threshold: ' num2str(PfaVector(minIdx))])
    disp('------------------------------------')
end
%     formatFig(gcf, [dataPath  'roc_CW_2'], 'en', figProp);
figure;
plot(JNRVector, minCminNMF)
hold on;
plot(JNRVector, minCmin);

%-------------------------------------------------------------------------
clear;
clc;

onset = 528;
offset = 4623;
window_length = 98;

thresholdVector = 0.1:0.05:0.9;
window_median_length_vector = 51:50:401;
monteCarloLoops = 100;
JNRVector = -20:0;

bandwidthVector = 10.72e6;
periodVector = 8.72e-6;

averageFP = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
averageTN = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
averageTP = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
averageFN = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));

averageTPR = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
averageFPR = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
stdTPR = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
stdFPR = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));

averageACC = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
stdACC = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
averageF1 = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
stdF1 = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
averagePrec = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
stdPrec = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
averageCmin = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));
stdCmin = zeros(length(JNRVector), length(thresholdVector), length(window_median_length_vector));

cminfun = @(tpr, fpr) sqrt(fpr.^2 + (1-tpr).^2);

for JNRIndex = 1:length(JNRVector)
    load(['results08_' num2str(JNRIndex)]);
    for bandwidthVectorIndex = 1:length(bandwidthVector)
        for periodIndex = 1:length(periodVector)
            for thresholdIndex = 1:length(thresholdVector)
                for window_median_length_index = 1:length(window_median_length_vector)
                    
                    fp = zeros(monteCarloLoops, 1);
                    tp = zeros(monteCarloLoops, 1);
                    fn = zeros(monteCarloLoops, 1);
                    tn = zeros(monteCarloLoops, 1);
                    
                    for loopIndex = 1:monteCarloLoops
                        fpRegion1 = squeeze(detection_res(loopIndex, bandwidthVectorIndex, periodIndex, 1, thresholdIndex, window_median_length_index, 1:onset));
                        tpRegion = squeeze(detection_res(loopIndex, bandwidthVectorIndex, periodIndex, 1, thresholdIndex, window_median_length_index, onset  + round(window_length*2):offset - round(window_length*2)));
                        fpRegion2 = squeeze(detection_res(loopIndex, bandwidthVectorIndex, periodIndex, 1, thresholdIndex, window_median_length_index, offset + round(window_length*2):end));
                        fpRegion = [fpRegion1;fpRegion2];
                        
                        fp(loopIndex) = sum(fpRegion);
                        tn(loopIndex) = length(fpRegion) - fp(loopIndex);
                        
                        tp(loopIndex) = sum(tpRegion);
                        fn(loopIndex) = length(tpRegion) - tp(loopIndex);
                    end
                    
                    acc = (tp + tn)./(tp + fn + tn + fp);
                    f1Score = 2*tp./(2*tp + fp + fn);
                    precision = tp./(tp+fp);
                    precision(isnan(precision)) = 0;
                    tpr = tp./(tp+fn);
                    fpr = fp./(fp+tn);
                    
                    cmin = cminfun(tpr, fpr);
                    
                    fpStd = std(fp);
                    tnStd = std(tn);
                    tpStd = std(tp);
                    fnStd = std(fn);
                    
                    averageFP(JNRIndex, thresholdIndex, window_median_length_index) = mean(fp);
                    averageTN(JNRIndex, thresholdIndex, window_median_length_index) = mean(tn);
                    averageTP(JNRIndex, thresholdIndex, window_median_length_index) = mean(tp);
                    averageFN(JNRIndex, thresholdIndex, window_median_length_index) = mean(fn);
                    
                    averageTPR(JNRIndex, thresholdIndex, window_median_length_index) = mean(tpr);
                    averageFPR(JNRIndex, thresholdIndex, window_median_length_index) = mean(fpr);
                    stdTPR(JNRIndex, thresholdIndex, window_median_length_index) = std(tpr);
                    stdFPR(JNRIndex, thresholdIndex, window_median_length_index) = std(fpr);
                    
                    averageACC(JNRIndex, thresholdIndex, window_median_length_index) = mean(acc);
                    stdACC(JNRIndex, thresholdIndex, window_median_length_index) = std(acc);
                    averageF1(JNRIndex, thresholdIndex, window_median_length_index) = mean(f1Score);
                    stdF1(JNRIndex, thresholdIndex, window_median_length_index) = std(f1Score);
                    averagePrec(JNRIndex, thresholdIndex, window_median_length_index) = mean(precision);
                    stdPrec(JNRIndex, thresholdIndex, window_median_length_index) = std(precision);
                    
                    averageCmin(JNRIndex, thresholdIndex, window_median_length_index) = mean(cmin);
                    stdCmin(JNRIndex, thresholdIndex, window_median_length_index) = std(cmin);
                end
            end
        end
    end
end

for i = 1:length(JNRVector)
    averageACCAux = squeeze(averageACC(i,:,:));

    [maxval, maxidx] = max(averageACCAux(:));
    [j, k] = ind2sub(size(averageACCAux), maxidx);
    disp('---------------------Accuracy-----------------------------------')
    disp('------------------------------------')
    disp(['-------------JNR: ' num2str(JNRVector(i)) ' dB-----------------------'])
    disp(['Accuracy: ' num2str(maxval) '+-' num2str(stdACC(i,j,k))])
    disp('------------------------------------')
    disp(['Pd: ' num2str(averageTPR(i,j,k)) '+-' num2str(stdTPR(i,j,k))])
    disp('------------------------------------')
    disp(['Pfa: ' num2str(averageFPR(i,j,k)) '+-' num2str(stdFPR(i,j,k))])
    disp('------------------------------------')
    disp(['Threshold: ' num2str(thresholdVector(j))])
    disp('------------------------------------')
    disp(['Window length: ' num2str(window_median_length_vector(k))])
    disp('------------------------------------')
    
%     plot(averageFPR(i,:,k), averageTPR(i,:,k))
    hold on    
end

for i = 1:length(JNRVector)
    averageF1Aux = squeeze(averageF1(i,:,:));

    [maxval, maxidx] = max(averageF1Aux(:));
    [j, k] = ind2sub(size(averageF1Aux), maxidx);
    disp('---------------------F1-score-----------------------------------')
    disp('------------------------------------')
    disp(['-------------JNR: ' num2str(JNRVector(i)) ' dB-----------------------'])
    disp(['F1-score: ' num2str(maxval) '+-' num2str(stdF1(i,j,k))])
    disp('------------------------------------')
    disp(['Pd: ' num2str(averageTPR(i,j,k)) '+-' num2str(stdTPR(i,j,k))])
    disp('------------------------------------')
    disp(['Pfa: ' num2str(averageFPR(i,j,k)) '+-' num2str(stdFPR(i,j,k))])
    disp('------------------------------------')
    disp(['Threshold: ' num2str(thresholdVector(j))])
    disp('------------------------------------')
    disp(['Window length: ' num2str(window_median_length_vector(k))])
    disp('------------------------------------')
    
%     plot(averageFPR(i,:,k), averageTPR(i,:,k))
end


for i = 1:length(JNRVector)
    averagePrecAux = squeeze(averagePrec(i,:,:));

    [maxval, maxidx] = max(averagePrecAux(:), [], 'includenan');
    [j, k] = ind2sub(size(averagePrecAux), maxidx);
    disp('---------------------Precision-----------------------------------')
    disp('------------------------------------')
    disp(['-------------JNR: ' num2str(JNRVector(i)) ' dB-----------------------'])
    disp(['Precision: ' num2str(maxval) '+-' num2str(stdPrec(i,j,k))])
    disp('------------------------------------')
    disp(['Pd: ' num2str(averageTPR(i,j,k)) '+-' num2str(stdTPR(i,j,k))])
    disp('------------------------------------')
    disp(['Pfa: ' num2str(averageFPR(i,j,k)) '+-' num2str(stdFPR(i,j,k))])
    disp('------------------------------------')
    disp(['Threshold: ' num2str(thresholdVector(j))])
    disp('------------------------------------')
    disp(['Window length: ' num2str(window_median_length_vector(k))])
    disp('------------------------------------')
    
%     plot(averageFPR(i,:,k), averageTPR(i,:,k))
end

for i = 1:length(JNRVector)

    averageCminAux = squeeze(averageCmin(i,:,:));
    [minCminNMF(i), minIdx] = min(averageCminAux(:));
    [j, k] = ind2sub(size(averageCminAux), minIdx);

    disp('---------------------Cmin-----------------------------------')
    disp('------------------------------------')
    disp(['-------------JNR: ' num2str(JNRVector(i)) ' dB-----------------------'])
    disp(['Cmin: ' num2str(minCminNMF(i)) '+-' num2str(stdCmin(i,j,k))])
    disp('------------------------------------')
    disp(['Pd: ' num2str(averageTPR(i,j,k)) '+-' num2str(stdTPR(i,j,k))])
    disp('------------------------------------')
    disp(['Pfa: ' num2str(averageFPR(i,j,k)) '+-' num2str(stdFPR(i,j,k))])
    disp('------------------------------------')
    disp(['Threshold: ' num2str(thresholdVector(j))])
    disp('------------------------------------')
    disp(['Window length: ' num2str(window_median_length_vector(k))])
    disp('------------------------------------')
end

clc;

addpath(['.' filesep 'data']);

load('resultsPai05.mat');

fs = 32.768e6;

numberOfRawSamples = 4096;
silenceSamples = round(20e-6*fs);
totalSamples = numberOfRawSamples + silenceSamples*2;
WinLBlock = 3;
MBlock = fix(totalSamples/WinLBlock);

JNRVector = -20:0;

PfaVector = logspace(-8, 0, 17);
frequencySliceCW = 400;
aux = 1:MBlock;
aux(frequencySliceCW) = 0;
aux = aux(aux~=0);

for JNRIndex = 1:length(JNRVector)
   
    for i = 1:length(PfaVector)
        x = squeeze(GoFBlockDeteflag(1, 1, JNRIndex, i,:,:));
        tpRegion = x(frequencySliceCW,:);
        fpRegion = x(aux,:);
        for j = 1:size(x, 2)
            tp(JNRIndex,j,i) = tpRegion(j);
            fn(JNRIndex,j,i) = 1 - tp(JNRIndex,j,i);
            fp(JNRIndex,j,i) = sum(fpRegion(:,j));
            tn(JNRIndex,j,i) = size(fpRegion(:,j), 1) - fp(JNRIndex,j,i);
        end
    end

end

tpr = tp./(tp+fn);
fpr = fp./(fp+tn);

tprAverage = squeeze(mean(tpr, 2));
fprAverage = squeeze(mean(fpr, 2));
tprStd = squeeze(std(tpr, [], 2));
fprStd = squeeze(std(fpr, [], 2));

cminfun = @(tpr, fpr) sqrt(fpr.^2 + (1-tpr).^2);
cmin = cminfun(tpr, fpr);
cminAverage = squeeze(mean(cmin, 2));
cminStd = squeeze(std(cmin, [], 2));    
    

for i = 1:length(JNRVector)
    [minCmin(i), minIdx] = min(cminAverage(i,:));

    disp('---------------------Cmin-----------------------------------')
    disp(['JNR: ' num2str(JNRVector(i)) ' dB'])
    disp('------------------------------------')
    disp(['Cmin: ' num2str(minCmin(i)) '+-' num2str(cminStd(i,minIdx))])
    disp('------------------------------------')
    disp(['Pd: ' num2str(tprAverage(i,minIdx)) '+-' num2str(tprStd(i,minIdx))])
    disp('------------------------------------')
    disp(['Pfa: ' num2str(fprAverage(i,minIdx)) '+-' num2str(fprStd(i,minIdx))])
    disp('------------------------------------')
    disp(['Threshold: ' num2str(PfaVector(minIdx))])
    disp('------------------------------------')
end
%     formatFig(gcf, [dataPath  'roc_CW_2'], 'en', figProp);
plot(JNRVector, minCminNMF)
hold on;
plot(JNRVector, minCmin);

legend('Dot, CW', 'Statistical, CW', 'Dot, Chirp', 'Statistical, Chirp');  


rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])
rmpath(['..' filesep '..' filesep '.' filesep 'data' filesep '08-12']);  
rmpath(['.' filesep 'data']);