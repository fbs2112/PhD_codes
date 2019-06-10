clear;
clc;
close all;

addpath(['..' filesep '.' filesep 'Misc'])
addpath(['..' filesep '.' filesep 'data' filesep '06-11']);    

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
dataPath = ['..' filesep '.' filesep 'figs' filesep '06-11' filesep];

threshold = 0.1:0.1:0.9;
stdVector = 7:2:13;


% load results01.mat;
% 
% figure;
% for i = 1:5
%     tpVector = squeeze(tp(:,i,1,:));
%     fpVector = squeeze(fp(:,i,1,:));
%     fnVector = squeeze(fn(:,i,1,:));
%     tnVector = squeeze(tn(:,i,1,:));
%    
%     accuracyStd = std((tpVector + tnVector)./(tpVector + fnVector + tnVector + fpVector), 1);
%     
%     tpr = tpVector./(tpVector+fnVector);
%     fpr = fpVector./(fpVector+tnVector);
%     
%     tprStd = std(tpr, 1);
%     fprStd = std(fpr, 1);
%     
%     tpMean = mean(tpVector, 1);
%     fpMean = mean(fpVector, 1);
%     fnMean = mean(fnVector, 1);
%     tnMean = mean(tnVector, 1);
%     
%     tpr = mean(tpr, 1);
%     fpr = mean(fpr, 1);
%     
%     plot(fpr, tpr);
%     hold on
%     
%     accuracy2(:,i) = (tpMean + tnMean)./(tpMean + fnMean + tnMean + fpMean);
%     [maxAccuracy(i), indexes(i)] = max(accuracy2(:,i), [], 1);
%     
%     disp('------------------------------------')
%     disp(['TPR: ' num2str(tpr(:, indexes(i))) '+' num2str(tprStd(:, indexes(i))) ' FPR: ' num2str(fpr(:, indexes(i))) '+' num2str(fprStd(:, indexes(i)))])
%     disp('------------------------------------')
%     disp(['Accuracy: ' num2str(accuracy2(indexes(i),i)) '+' num2str(accuracyStd(indexes(i)))])
%     disp('------------------------------------')
%     disp(['Threshold: ' num2str(threshold(indexes(i)))])
%     disp('------------------------------------')
%     
% end
% 
% xAux = linspace(0, 1, numel(fpr));
% yAux = linspace(0, 1, numel(tpr));
% plot(xAux, yAux, '--');
% box 'off'
% grid on
% xlabel('Probability of false alarm');
% ylabel('Probability of detection');
% legend('-20 dB', '-15 dB', '-10 dB', '-5 dB', '0 dB', 'Random guess', 'Location', 'southeast');
% 
% expName = 'detect_2chirp_sim';
% 
% save_fig = false;
% if save_fig
%     formatFig(gcf, [dataPath expName '_' 'roc_'  '1_'], 'en', figProp); %#ok<*UNRCH>
% end
% %-------------------------------------------------------
% 
% clc;
% 
% load results02.mat;
% 
% figure;
% for i = 1:5
%     tpVector = squeeze(tp(:,i,1,:));
%     fpVector = squeeze(fp(:,i,1,:));
%     fnVector = squeeze(fn(:,i,1,:));
%     tnVector = squeeze(tn(:,i,1,:));
%    
%     accuracyStd = std((tpVector + tnVector)./(tpVector + fnVector + tnVector + fpVector), 1);
%     
%     tpr = tpVector./(tpVector+fnVector);
%     fpr = fpVector./(fpVector+tnVector);
%     
%     tprStd = std(tpr, 1);
%     fprStd = std(fpr, 1);
%     
%     tpMean = mean(tpVector, 1);
%     fpMean = mean(fpVector, 1);
%     fnMean = mean(fnVector, 1);
%     tnMean = mean(tnVector, 1);
%     
%     tpr = mean(tpr, 1);
%     fpr = mean(fpr, 1);
%     
%     plot(fpr, tpr);
%     hold on
%     
%     accuracy2(:,i) = (tpMean + tnMean)./(tpMean + fnMean + tnMean + fpMean);
%     [maxAccuracy(i), indexes(i)] = max(accuracy2(:,i), [], 1);
%     
%     disp('------------------------------------')
%     disp(['TPR: ' num2str(tpr(:, indexes(i))) '+' num2str(tprStd(:, indexes(i))) ' FPR: ' num2str(fpr(:, indexes(i))) '+' num2str(fprStd(:, indexes(i)))])
%     disp('------------------------------------')
%     disp(['Accuracy: ' num2str(accuracy2(indexes(i),i)) '+' num2str(accuracyStd(indexes(i)))])
%     disp('------------------------------------')
%     disp(['Threshold: ' num2str(threshold(indexes(i)))])
%     disp('------------------------------------')
%     
% end
% 
% xAux = linspace(0, 1, numel(fpr));
% yAux = linspace(0, 1, numel(tpr));
% plot(xAux, yAux, '--');
% box 'off'
% grid on
% xlabel('Probability of false alarm');
% ylabel('Probability of detection');
% legend('-20 dB', '-15 dB', '-10 dB', '-5 dB', '0 dB', 'Random guess', 'Location', 'southeast');
% 
% expName = 'detect_2chirp_sim';
% 
% save_fig = false;
% if save_fig
%     formatFig(gcf, [dataPath expName '_' 'roc_'  '1_'], 'en', figProp); %#ok<*UNRCH>
% end
% %-------------------------------------------------------
% 
% clc;

% load resultsEdiz.mat;
% expName = 'output_analysis';

JNR = [-20 -15 -10 -5 0 10];

% for JNRIndex = 1:length(JNR)
% 
%     figure;
%     plot(t*1e6, output(:, JNRIndex, 1)./ max(output(:, JNRIndex, 1)))
%     ylabel('Normalized Magnitude');
%     xlabel('Time [$\mu$s]');
%     ylim([0 1.1])
%     xlim([min(t) max(t)]*1e6);
%     formatFig(gcf, [dataPath expName '_' 'sim_' num2str(JNR(JNRIndex))], 'en', figProp);
% 
%     figure
%     plot(t*1e6, outputVar(:, JNRIndex, 1) ./ max(outputVar(:, JNRIndex, 1)))
%     ylabel('Normalized Magnitude');
%     xlabel('Time [$\mu$s]');
%     ylim([0 1.1])
%     xlim([min(t) max(t)]*1e6);
%     formatFig(gcf, [dataPath expName '_' 'sim_var_' num2str(JNR(JNRIndex))], 'en', figProp);
% 
%     figure
%     plot(t*1e6, outputTK(:, JNRIndex, 1))
%     ylabel('Normalized Magnitude');
%     xlabel('Time [$\mu$s]');
%     % ylim([0 1.1])
%     xlim([min(t) max(t)]*1e6);
%     formatFig(gcf, [dataPath expName '_' 'TK_' num2str(JNR(JNRIndex))], 'en', figProp);
% 
%     figure
%     plot(t*1e6, outputTKVar(:, JNRIndex, 1))
%     ylabel('Magnitude');
%     xlabel('Time [$\mu$s]');
%     % ylim([0 1.1])
%     xlim([min(t) max(t)]*1e6);
%     formatFig(gcf, [dataPath expName '_' 'TK_var_' num2str(JNR(JNRIndex))], 'en', figProp);
% end

close all;

load resultsEdiz2.mat;
expName = 'output_analysis_signal';
onset = 528;
offset = 3911;

for JNRIndex = 1:length(JNR)

    figure;
    plot(t*1e6, output(:, JNRIndex, 1)./ max(output(:, JNRIndex, 1)))
    ylabel('Normalized Magnitude');
    xlabel('Time [$\mu$s]');
    
    hold on;
    line([t(onset) t(onset)]*1e6, [0 1.1], 'Color','black','LineStyle','--');
    line([t(offset) t(offset)]*1e6, [0 1.1], 'Color','black','LineStyle','--');
    ylim([0 1.1])
    xlim([min(t) max(t)]*1e6);
    formatFig(gcf, [dataPath expName '_' 'sim_' num2str(JNR(JNRIndex))], 'en', figProp);

    figure
    plot(t*1e6, outputVar(:, JNRIndex, 1) ./ max(outputVar(:, JNRIndex, 1)))
    ylabel('Normalized Magnitude');
    xlabel('Time [$\mu$s]');
    hold on;
    line([t(onset) t(onset)]*1e6, [0 1.1], 'Color','black','LineStyle','--');
    line([t(offset) t(offset)]*1e6, [0 1.1], 'Color','black','LineStyle','--');
    ylim([0 1.1])
    xlim([min(t) max(t)]*1e6);
    formatFig(gcf, [dataPath expName '_' 'sim_var_' num2str(JNR(JNRIndex))], 'en', figProp);

    figure
    plot(t*1e6, outputTK(:, JNRIndex, 1))
    ylabel('Magnitude');
    xlabel('Time [$\mu$s]');
    hold on;
    line([t(onset) t(onset)]*1e6, [1.1*min(outputTK(:, JNRIndex, 1)) 1.1*max(outputTK(:, JNRIndex, 1))], 'Color','black','LineStyle','--');
    line([t(offset) t(offset)]*1e6, [1.1*min(outputTK(:, JNRIndex, 1)) 1.1*max(outputTK(:, JNRIndex, 1))], 'Color','black','LineStyle','--');
    ylim([1.1*min(outputTK(:, JNRIndex, 1)) 1.1*max(outputTK(:, JNRIndex, 1))])
    xlim([min(t) max(t)]*1e6);
    formatFig(gcf, [dataPath expName '_' 'TK_' num2str(JNR(JNRIndex))], 'en', figProp);

    figure
    plot(t*1e6, outputTKVar(:, JNRIndex, 1))
    ylabel('Magnitude');
    xlabel('Time [$\mu$s]');
    hold on;
    line([t(onset) t(onset)]*1e6, [1.1*min(outputTKVar(:, JNRIndex, 1)) 1.1*max(outputTKVar(:, JNRIndex, 1))], 'Color','black','LineStyle','--');
    line([t(offset) t(offset)]*1e6, [1.1*min(outputTKVar(:, JNRIndex, 1)) 1.1*max(outputTKVar(:, JNRIndex, 1))], 'Color','black','LineStyle','--');
    ylim([1.1*min(outputTKVar(:, JNRIndex, 1)) 1.1*max(outputTKVar(:, JNRIndex, 1))])
    xlim([min(t) max(t)]*1e6);
    formatFig(gcf, [dataPath expName '_' 'TK_var_' num2str(JNR(JNRIndex))], 'en', figProp);
end


rmpath(['..' filesep '.' filesep 'data' filesep '06-11']);    
rmpath(['..' filesep '.' filesep 'Misc'])
