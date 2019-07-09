clear;
clc;
close all;

addpath(['..' filesep '.' filesep 'Misc'])
addpath(['..' filesep '.' filesep 'data' filesep '07-09']);    

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
dataPath = ['..' filesep '.' filesep 'figs' filesep '06-11' filesep];

lengthVector = 200:100:500;
stdVector2 = [2 3 4];
forgettingFactorVector = [0.3 0.4 0.5 0.6];
monteCarloLoops = 100;
JNRVector = [-10 -5 0];
onset = 528;
offset = 3911 + 100;

numberOfTruePositives = 13;

averageFP = zeros(length(JNRVector), length(lengthVector), length(stdVector2), length(forgettingFactorVector));
averageTN = zeros(length(JNRVector), length(lengthVector), length(stdVector2), length(forgettingFactorVector));
averageTP = zeros(length(JNRVector), length(lengthVector), length(stdVector2), length(forgettingFactorVector));
averageFN = zeros(length(JNRVector), length(lengthVector), length(stdVector2), length(forgettingFactorVector));

averageTPR = zeros(length(JNRVector), length(lengthVector), length(stdVector2), length(forgettingFactorVector));
averageFPR = zeros(length(JNRVector), length(lengthVector), length(stdVector2), length(forgettingFactorVector));
stdTPR = zeros(length(JNRVector), length(lengthVector), length(stdVector2), length(forgettingFactorVector));
stdFPR = zeros(length(JNRVector), length(lengthVector), length(stdVector2), length(forgettingFactorVector));

averageACC = zeros(length(JNRVector), length(lengthVector), length(stdVector2), length(forgettingFactorVector));
stdACC = zeros(length(JNRVector), length(lengthVector), length(stdVector2), length(forgettingFactorVector));
averageF1 = zeros(length(JNRVector), length(lengthVector), length(stdVector2), length(forgettingFactorVector));
stdF1 = zeros(length(JNRVector), length(lengthVector), length(stdVector2), length(forgettingFactorVector));
averagePrec = zeros(length(JNRVector), length(lengthVector), length(stdVector2), length(forgettingFactorVector));
stdPrec = zeros(length(JNRVector), length(lengthVector), length(stdVector2), length(forgettingFactorVector));


load results01.mat;

for JNRIndex = 1:length(JNRVector)
    for lengthVectorIndex = 1:length(lengthVector)
        for stdVector2Index = 1:length(stdVector2)
            for forgettingFactorVectorIndex = 1:length(forgettingFactorVector)
                
                fpRegion1 = squeeze(detection_res(:, JNRIndex, 1, lengthVectorIndex, stdVector2Index, forgettingFactorVectorIndex, 1:onset));
                fpRegion2 = squeeze(detection_res(:, JNRIndex, 1, lengthVectorIndex, stdVector2Index, forgettingFactorVectorIndex, offset:end));
                tpRegion = squeeze(detection_res(:, JNRIndex, 1, lengthVectorIndex, stdVector2Index, forgettingFactorVectorIndex, onset+1:offset-1));
                
                fp = zeros(monteCarloLoops, 1);
                tp = zeros(monteCarloLoops, 1);
                fn = zeros(monteCarloLoops, 1);
                tn = zeros(monteCarloLoops, 1);
                
                for loopIndex = 1:monteCarloLoops
                    
                    numberOfPeaks1 = length(findpeaks(fpRegion1(loopIndex,:)));
                    numberOfPeaks2 = length(findpeaks(fpRegion2(loopIndex,:)));
                    numberOfPeaks3 = length(findpeaks(tpRegion(loopIndex,:)));
                    
                    if numberOfPeaks1 > 1
                        fp(loopIndex) = fp(loopIndex) + 1;
                    else
                        tn(loopIndex) = tn(loopIndex) + 1;
                    end
                    
                    if numberOfPeaks2 > 1
                        fp(loopIndex) = fp(loopIndex) + 1;
                    else
                        tn(loopIndex) = tn(loopIndex) + 1;
                    end
                    
                    if numberOfPeaks3 == numberOfTruePositives
                        tp(loopIndex) = numberOfTruePositives;
                    end
                        
                    if numberOfPeaks3 > numberOfTruePositives
                        tp(loopIndex) = numberOfTruePositives;
                        fp(loopIndex) = numberOfPeaks3 - numberOfTruePositives;
                    end
                    
                    if numberOfPeaks3 < numberOfTruePositives
                        tp(loopIndex) = numberOfPeaks3;
                        fn(loopIndex) = numberOfTruePositives - numberOfPeaks3;
                    end
                    
                end
                
                acc = (tp + tn)./(tp + fn + tn + fp);
                f1Score = 2*tp./(2*tp + fp + fn);
                precision = tp./(tp+fp);
                tpr = tp./(tp+fn);
                fpr = fp./(fp+tn);
                
                fpStd = std(fp);
                tnStd = std(tn);
                tpStd = std(tp);
                fnStd = std(fn);
                
                averageFP(JNRIndex, lengthVectorIndex, stdVector2Index, forgettingFactorVectorIndex) = mean(fp);
                averageTN(JNRIndex, lengthVectorIndex, stdVector2Index, forgettingFactorVectorIndex) = mean(tn);
                averageTP(JNRIndex, lengthVectorIndex, stdVector2Index, forgettingFactorVectorIndex) = mean(tp);
                averageFN(JNRIndex, lengthVectorIndex, stdVector2Index, forgettingFactorVectorIndex) = mean(fn);
                
                averageTPR(JNRIndex, lengthVectorIndex, stdVector2Index, forgettingFactorVectorIndex) = mean(tpr);
                averageFPR(JNRIndex, lengthVectorIndex, stdVector2Index, forgettingFactorVectorIndex) = mean(fpr);
                stdTPR(JNRIndex, lengthVectorIndex, stdVector2Index, forgettingFactorVectorIndex) = std(tpr);
                stdFPR(JNRIndex, lengthVectorIndex, stdVector2Index, forgettingFactorVectorIndex) = std(fpr);
                
                averageACC(JNRIndex, lengthVectorIndex, stdVector2Index, forgettingFactorVectorIndex) = mean(acc);
                stdACC(JNRIndex, lengthVectorIndex, stdVector2Index, forgettingFactorVectorIndex) = std(acc);
                averageF1(JNRIndex, lengthVectorIndex, stdVector2Index, forgettingFactorVectorIndex) = mean(f1Score);
                stdF1(JNRIndex, lengthVectorIndex, stdVector2Index, forgettingFactorVectorIndex) = std(f1Score);
                averagePrec(JNRIndex, lengthVectorIndex, stdVector2Index, forgettingFactorVectorIndex) = mean(precision);
                stdPrec(JNRIndex, lengthVectorIndex, stdVector2Index, forgettingFactorVectorIndex) = std(precision);
            end
        end
    end
end


for i = 1:1%length(JNRVector)
    averageACCAux = squeeze(averageACC(i,:,:,:));

    [maxval, maxidx] = max(averageACCAux(:));
    [j, k, l] = ind2sub(size(averageACCAux), maxidx);

    disp(['-------------JNR: ' num2str(JNRVector(i)) ' dB-----------------------'])
    disp(['Accuracy: ' num2str(maxval) '+-' num2str(stdACC(i,j,k,l))])
    disp('------------------------------------')
    disp(['Pd: ' num2str(averageTPR(i,j,k,l)) '+-' num2str(stdTPR(i,j,k,l))])
    disp('------------------------------------')
    disp(['Pfa: ' num2str(averageFPR(i,j,k,l)) '+-' num2str(stdFPR(i,j,k,l))])
    disp('------------------------------------')
    disp(['Length: ' num2str(lengthVector(j))])
    disp('------------------------------------')
    disp(['Std: ' num2str(stdVector2(k))])
    disp('------------------------------------')
    disp(['Forgetting factor: ' num2str(forgettingFactorVector(l))])
    disp('------------------------------------')
end


for i = 1:1%length(JNRVector)
    averageF1Aux = squeeze(averageF1(i,:,:,:));

    [maxval, maxidx] = max(averageF1Aux(:));
    [j, k, l] = ind2sub(size(averageF1Aux), maxidx);

    disp(['-------------JNR: ' num2str(JNRVector(i)) ' dB-----------------------'])
    disp(['F1 score: ' num2str(maxval) '+-' num2str(stdF1(i,j,k,l))])
    disp('------------------------------------')
    disp(['Pd: ' num2str(averageTPR(i,j,k,l)) '+-' num2str(stdTPR(i,j,k,l))])
    disp('------------------------------------')
    disp(['Pfa: ' num2str(averageFPR(i,j,k,l)) '+-' num2str(stdFPR(i,j,k,l))])
    disp('------------------------------------')
    disp(['Length: ' num2str(lengthVector(j))])
    disp('------------------------------------')
    disp(['Std: ' num2str(stdVector2(k))])
    disp('------------------------------------')
    disp(['Forgetting factor: ' num2str(forgettingFactorVector(l))])
    disp('------------------------------------')
end


for i = 1:1%length(JNRVector)
    averagePrecAux = squeeze(averagePrec(i,:,:,:));

    [maxval, maxidx] = max(averagePrecAux(:));
    [j, k, l] = ind2sub(size(averagePrecAux), maxidx);

    disp(['-------------JNR: ' num2str(JNRVector(i)) ' dB-----------------------'])
    disp(['Precision: ' num2str(maxval) '+-' num2str(stdPrec(i,j,k,l))])
    disp('------------------------------------')
    disp(['Pd: ' num2str(averageTPR(i,j,k,l)) '+-' num2str(stdTPR(i,j,k,l))])
    disp('------------------------------------')
    disp(['Pfa: ' num2str(averageFPR(i,j,k,l)) '+-' num2str(stdFPR(i,j,k,l))])
    disp('------------------------------------')
    disp(['Length: ' num2str(lengthVector(j))])
    disp('------------------------------------')
    disp(['Std: ' num2str(stdVector2(k))])
    disp('------------------------------------')
    disp(['Forgetting factor: ' num2str(forgettingFactorVector(l))])
    disp('------------------------------------')
end


rmpath(['..' filesep '.' filesep 'data' filesep '07-09']);    
rmpath(['..' filesep '.' filesep 'Misc'])

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



