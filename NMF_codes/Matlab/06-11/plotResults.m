clear;
clc;
close all;

addpath(['..' filesep '.' filesep 'Misc'])


linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
dataPath = ['..' filesep '.' filesep 'figs' filesep '05-23' filesep];

threshold = 0.1:0.1:0.9;
stdVector = 7:2:13;

addpath(['..' filesep '.' filesep 'data' filesep '06-11']);    

load results01.mat;

figure;
for i = 1:5
    tpVector = squeeze(tp(:,i,1,:));
    fpVector = squeeze(fp(:,i,1,:));
    fnVector = squeeze(fn(:,i,1,:));
    tnVector = squeeze(tn(:,i,1,:));
   
    accuracyStd = std((tpVector + tnVector)./(tpVector + fnVector + tnVector + fpVector), 1);
    
    tpr = tpVector./(tpVector+fnVector);
    fpr = fpVector./(fpVector+tnVector);
    
    tprStd = std(tpr, 1);
    fprStd = std(fpr, 1);
    
    tpMean = mean(tpVector, 1);
    fpMean = mean(fpVector, 1);
    fnMean = mean(fnVector, 1);
    tnMean = mean(tnVector, 1);
    
    tpr = mean(tpr, 1);
    fpr = mean(fpr, 1);
    
    plot(fpr, tpr);
    hold on
    
    accuracy2(:,i) = (tpMean + tnMean)./(tpMean + fnMean + tnMean + fpMean);
    [maxAccuracy(i), indexes(i)] = max(accuracy2(:,i), [], 1);
    
    disp('------------------------------------')
    disp(['TPR: ' num2str(tpr(:, indexes(i))) '+' num2str(tprStd(:, indexes(i))) ' FPR: ' num2str(fpr(:, indexes(i))) '+' num2str(fprStd(:, indexes(i)))])
    disp('------------------------------------')
    disp(['Accuracy: ' num2str(accuracy2(indexes(i),i)) '+' num2str(accuracyStd(indexes(i)))])
    disp('------------------------------------')
    disp(['Threshold: ' num2str(threshold(indexes(i)))])
    disp('------------------------------------')
    
end

xAux = linspace(0, 1, numel(fpr));
yAux = linspace(0, 1, numel(tpr));
plot(xAux, yAux, '--');
box 'off'
grid on
xlabel('Probability of false alarm');
ylabel('Probability of detection');
legend('-20 dB', '-15 dB', '-10 dB', '-5 dB', '0 dB', 'Random guess', 'Location', 'southeast');

expName = 'detect_2chirp_sim';

save_fig = false;
if save_fig
    formatFig(gcf, [dataPath expName '_' 'roc_'  '1_'], 'en', figProp); %#ok<*UNRCH>
end
%-------------------------------------------------------

clc;

load results02.mat;

figure;
for i = 1:5
    tpVector = squeeze(tp(:,i,1,:));
    fpVector = squeeze(fp(:,i,1,:));
    fnVector = squeeze(fn(:,i,1,:));
    tnVector = squeeze(tn(:,i,1,:));
   
    accuracyStd = std((tpVector + tnVector)./(tpVector + fnVector + tnVector + fpVector), 1);
    
    tpr = tpVector./(tpVector+fnVector);
    fpr = fpVector./(fpVector+tnVector);
    
    tprStd = std(tpr, 1);
    fprStd = std(fpr, 1);
    
    tpMean = mean(tpVector, 1);
    fpMean = mean(fpVector, 1);
    fnMean = mean(fnVector, 1);
    tnMean = mean(tnVector, 1);
    
    tpr = mean(tpr, 1);
    fpr = mean(fpr, 1);
    
    plot(fpr, tpr);
    hold on
    
    accuracy2(:,i) = (tpMean + tnMean)./(tpMean + fnMean + tnMean + fpMean);
    [maxAccuracy(i), indexes(i)] = max(accuracy2(:,i), [], 1);
    
    disp('------------------------------------')
    disp(['TPR: ' num2str(tpr(:, indexes(i))) '+' num2str(tprStd(:, indexes(i))) ' FPR: ' num2str(fpr(:, indexes(i))) '+' num2str(fprStd(:, indexes(i)))])
    disp('------------------------------------')
    disp(['Accuracy: ' num2str(accuracy2(indexes(i),i)) '+' num2str(accuracyStd(indexes(i)))])
    disp('------------------------------------')
    disp(['Threshold: ' num2str(threshold(indexes(i)))])
    disp('------------------------------------')
    
end

xAux = linspace(0, 1, numel(fpr));
yAux = linspace(0, 1, numel(tpr));
plot(xAux, yAux, '--');
box 'off'
grid on
xlabel('Probability of false alarm');
ylabel('Probability of detection');
legend('-20 dB', '-15 dB', '-10 dB', '-5 dB', '0 dB', 'Random guess', 'Location', 'southeast');

expName = 'detect_2chirp_sim';

save_fig = false;
if save_fig
    formatFig(gcf, [dataPath expName '_' 'roc_'  '1_'], 'en', figProp); %#ok<*UNRCH>
end
%-------------------------------------------------------

clc;

load resultsEdiz2.mat;

JNR = [-20 -15 -10 -5 0 10];
JNRIndex = 3;   

figure;

plot(output(:, JNRIndex, 1)./ max(output(:, JNRIndex, 1)))
figure
plot(outputVar(:, JNRIndex, 1) ./ max(outputVar(:, JNRIndex, 1)))
figure
plot(outputTK(:, JNRIndex, 1) ./ max(outputTK(:, JNRIndex, 1)))

figure

% x = filter(ones(100, 1), 1, outputTK(:, JNRIndex, 1) ./ max(outputTK(:, JNRIndex, 1)))

plot(envelope(outputTK(:, JNRIndex, 1) ./ max(outputTK(:, JNRIndex, 1))))


rmpath(['..' filesep '.' filesep 'data' filesep '06-11']);    
rmpath(['..' filesep '.' filesep 'Misc'])
