clear;
clc;
close all;

addpath(['..' filesep '.' filesep 'Misc'])
addpath(['..' filesep '.' filesep 'data' filesep '07-24']);    

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
dataPath = ['..' filesep '.' filesep 'figs' filesep '06-11' filesep];

lengthVector = 500;
stdVector2 = [2];
forgettingFactorVector = [0.3];
JNRVector = [0];

numberOfTruePositives = 13;

averageFP = zeros(length(JNRVector), length(lengthVector), length(stdVector2), length(forgettingFactorVector));
averageTN = zeros(length(JNRVector), length(lengthVector), length(stdVector2), length(forgettingFactorVector));

averageFPR = zeros(length(JNRVector), length(lengthVector), length(stdVector2), length(forgettingFactorVector));
stdFPR = zeros(length(JNRVector), length(lengthVector), length(stdVector2), length(forgettingFactorVector));


load results01.mat;

for JNRIndex = 1:length(JNRVector)
    for lengthVectorIndex = 1:length(lengthVector)
        for stdVector2Index = 1:length(stdVector2)
            for forgettingFactorVectorIndex = 1:length(forgettingFactorVector)
                
                fpRegion = squeeze(detection_res(:, JNRIndex, 1, lengthVectorIndex, stdVector2Index, forgettingFactorVectorIndex, 1:end));
                
                fp = sum(fpRegion, 2);
                tn = 3969 - fp;
                
                fpr = fp./(fp+tn);
                fpStd = std(fp);
                tnStd = std(tn);
                
                averageFP(JNRIndex, lengthVectorIndex, stdVector2Index, forgettingFactorVectorIndex) = mean(fp);
                averageTN(JNRIndex, lengthVectorIndex, stdVector2Index, forgettingFactorVectorIndex) = mean(tn);
                averageFPR(JNRIndex, lengthVectorIndex, stdVector2Index, forgettingFactorVectorIndex) = mean(fpr);
                stdFPR(JNRIndex, lengthVectorIndex, stdVector2Index, forgettingFactorVectorIndex) = std(fpr);
                
            end
        end
    end
end


for i = 1:1%length(JNRVector)

    j = 1;
    k = 1;
    l = 1;
    disp(['-------------JNR: ' num2str(JNRVector(i)) ' dB-----------------------'])
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
clc;

load results02.mat;

onset = 528;
offset = 3911;
window_length = 98;

thresholdVector = 0.1:0.05:0.9;
window_median_length_vector = 51:50:401;
monteCarloLoops = 100;
JNRVector = [-20 -15 -10 -5 0];

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

for JNRIndex = 1:length(JNRVector)
    for thresholdIndex = 1:length(thresholdVector)
        for window_median_length_index = 1:length(window_median_length_vector)
            
            fp = zeros(monteCarloLoops, 1);
            tp = zeros(monteCarloLoops, 1);
            fn = zeros(monteCarloLoops, 1);
            tn = zeros(monteCarloLoops, 1);
            
            for loopIndex = 1:monteCarloLoops
                fpRegion1 = squeeze(detection_res(loopIndex, JNRIndex, thresholdIndex, window_median_length_index, 1:onset));
                tpRegion = squeeze(detection_res(loopIndex, JNRIndex, thresholdIndex, window_median_length_index, onset  + round(window_length*2):offset - round(window_length*2)));
                fpRegion2 = squeeze(detection_res(loopIndex, JNRIndex, thresholdIndex, window_median_length_index, offset + round(window_length*2):end));
                
                if any(fpRegion1)
                    fp(loopIndex) = fp(loopIndex) + 1;
                else
                    tn(loopIndex) = tn(loopIndex) + 1;
                end
                
                if any(fpRegion2)
                    fp(loopIndex) = fp(loopIndex) + 1;
                else
                    tn(loopIndex) = tn(loopIndex) + 1;
                end
                
                if all(tpRegion)
                    tp(loopIndex) = tp(loopIndex) + 1;
                else
                    fn(loopIndex) = fn(loopIndex) + 1;
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

xAux = linspace(0, 1, numel(averageTPR(i, :, k)));
yAux = linspace(0, 1, numel(averageFPR(i, :, k)));
plot(xAux, yAux, '--');
box 'off'
grid on
xlabel('Probability of false alarm');
ylabel('Probability of detection');
legend('-20 dB', '-15 dB', '-10 dB', '-5 dB', '0 dB', 'Random guess', 'Location', 'southeast');


figure;
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
    hold on    
end

xAux = linspace(0, 1, numel(averageTPR(i, :, k)));
yAux = linspace(0, 1, numel(averageFPR(i, :, k)));
plot(xAux, yAux, '--');
box 'off'
grid on
xlabel('Probability of false alarm');
ylabel('Probability of detection');
legend('-20 dB', '-15 dB', '-10 dB', '-5 dB', '0 dB', 'Random guess', 'Location', 'southeast');


figure;
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
    hold on 
end

xAux = linspace(0, 1, numel(averageTPR(i, :, k)));
yAux = linspace(0, 1, numel(averageFPR(i, :, k)));
plot(xAux, yAux, '--');
box 'off'
grid on
xlabel('Probability of false alarm');
ylabel('Probability of detection');
legend('-20 dB', '-15 dB', '-10 dB', '-5 dB', '0 dB', 'Random guess', 'Location', 'southeast');


clc;

load results03.mat;

onset = 528;
offset = 3911;
window_length = 98;

thresholdVector = 0.1:0.05:0.9;
window_median_length_vector = 51:50:401;
monteCarloLoops = 100;
JNRVector = [-20 -15 -10 -5 0];

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

for JNRIndex = 1:length(JNRVector)
    for thresholdIndex = 1:length(thresholdVector)
        for window_median_length_index = 1:length(window_median_length_vector)
            
            fp = zeros(monteCarloLoops, 1);
            tp = zeros(monteCarloLoops, 1);
            fn = zeros(monteCarloLoops, 1);
            tn = zeros(monteCarloLoops, 1);
            
            for loopIndex = 1:monteCarloLoops
                fpRegion1 = squeeze(detection_res(loopIndex, JNRIndex, thresholdIndex, window_median_length_index, 1:onset));
                tpRegion = squeeze(detection_res(loopIndex, JNRIndex, thresholdIndex, window_median_length_index, onset  + round(window_length*2):offset - round(window_length*2)));
                fpRegion2 = squeeze(detection_res(loopIndex, JNRIndex, thresholdIndex, window_median_length_index, offset + round(window_length*2):end));
                
                if any(fpRegion1)
                    fp(loopIndex) = fp(loopIndex) + 1;
                else
                    tn(loopIndex) = tn(loopIndex) + 1;
                end
                
                if any(fpRegion2)
                    fp(loopIndex) = fp(loopIndex) + 1;
                else
                    tn(loopIndex) = tn(loopIndex) + 1;
                end
                
                if all(tpRegion)
                    tp(loopIndex) = tp(loopIndex) + 1;
                else
                    fn(loopIndex) = fn(loopIndex) + 1;
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

xAux = linspace(0, 1, numel(averageTPR(i, :, k)));
yAux = linspace(0, 1, numel(averageFPR(i, :, k)));
plot(xAux, yAux, '--');
box 'off'
grid on
xlabel('Probability of false alarm');
ylabel('Probability of detection');
legend('-20 dB', '-15 dB', '-10 dB', '-5 dB', '0 dB', 'Random guess', 'Location', 'southeast');


figure;
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
    hold on    
end

xAux = linspace(0, 1, numel(averageTPR(i, :, k)));
yAux = linspace(0, 1, numel(averageFPR(i, :, k)));
plot(xAux, yAux, '--');
box 'off'
grid on
xlabel('Probability of false alarm');
ylabel('Probability of detection');
legend('-20 dB', '-15 dB', '-10 dB', '-5 dB', '0 dB', 'Random guess', 'Location', 'southeast');


figure;
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
    hold on 
end

xAux = linspace(0, 1, numel(averageTPR(i, :, k)));
yAux = linspace(0, 1, numel(averageFPR(i, :, k)));
plot(xAux, yAux, '--');
box 'off'
grid on
xlabel('Probability of false alarm');
ylabel('Probability of detection');
legend('-20 dB', '-15 dB', '-10 dB', '-5 dB', '0 dB', 'Random guess', 'Location', 'southeast');

rmpath(['..' filesep '.' filesep 'data' filesep '07-24']);    
rmpath(['..' filesep '.' filesep 'Misc'])

%%
%Count true positives, negatives, false positives, negatives sample by
%sample
    
clear;
clc;
close all;

addpath(['..' filesep '.' filesep 'Misc'])
addpath(['..' filesep '.' filesep 'data' filesep '07-24']);    

load results02.mat;

onset = 528;
offset = 3911;
window_length = 98;

thresholdVector = 0.1:0.05:0.9;
window_median_length_vector = 51:50:401;
monteCarloLoops = 100;
JNRVector = [-20 -15 -10 -5 0];

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

for JNRIndex = 1:length(JNRVector)
    for thresholdIndex = 1:length(thresholdVector)
        for window_median_length_index = 1:length(window_median_length_vector)
            
            fp = zeros(monteCarloLoops, 1);
            tp = zeros(monteCarloLoops, 1);
            fn = zeros(monteCarloLoops, 1);
            tn = zeros(monteCarloLoops, 1);
            
            for loopIndex = 1:monteCarloLoops
                fpRegion1 = squeeze(detection_res(loopIndex, JNRIndex, thresholdIndex, window_median_length_index, 1:onset));
                tpRegion = squeeze(detection_res(loopIndex, JNRIndex, thresholdIndex, window_median_length_index, onset  + round(window_length*2):offset - round(window_length*2)));
                fpRegion2 = squeeze(detection_res(loopIndex, JNRIndex, thresholdIndex, window_median_length_index, offset + round(window_length*2):end));
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

xAux = linspace(0, 1, numel(averageTPR(i, :, k)));
yAux = linspace(0, 1, numel(averageFPR(i, :, k)));
plot(xAux, yAux, '--');
box 'off'
grid on
xlabel('Probability of false alarm');
ylabel('Probability of detection');
legend('-20 dB', '-15 dB', '-10 dB', '-5 dB', '0 dB', 'Random guess', 'Location', 'southeast');


figure;
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
    hold on    
end

xAux = linspace(0, 1, numel(averageTPR(i, :, k)));
yAux = linspace(0, 1, numel(averageFPR(i, :, k)));
plot(xAux, yAux, '--');
box 'off'
grid on
xlabel('Probability of false alarm');
ylabel('Probability of detection');
legend('-20 dB', '-15 dB', '-10 dB', '-5 dB', '0 dB', 'Random guess', 'Location', 'southeast');


figure;
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
    hold on 
end

xAux = linspace(0, 1, numel(averageTPR(i, :, k)));
yAux = linspace(0, 1, numel(averageFPR(i, :, k)));
plot(xAux, yAux, '--');
box 'off'
grid on
xlabel('Probability of false alarm');
ylabel('Probability of detection');
legend('-20 dB', '-15 dB', '-10 dB', '-5 dB', '0 dB', 'Random guess', 'Location', 'southeast');



load results03.mat;

onset = 528;
offset = 3911;
window_length = 98;

thresholdVector = 0.1:0.05:0.9;
window_median_length_vector = 51:50:401;
monteCarloLoops = 100;
JNRVector = [-20 -15 -10 -5 0];

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

for JNRIndex = 1:length(JNRVector)
    for thresholdIndex = 1:length(thresholdVector)
        for window_median_length_index = 1:length(window_median_length_vector)
            
            fp = zeros(monteCarloLoops, 1);
            tp = zeros(monteCarloLoops, 1);
            fn = zeros(monteCarloLoops, 1);
            tn = zeros(monteCarloLoops, 1);
            
            for loopIndex = 1:monteCarloLoops
                fpRegion1 = squeeze(detection_res(loopIndex, JNRIndex, thresholdIndex, window_median_length_index, 1:onset));
                tpRegion = squeeze(detection_res(loopIndex, JNRIndex, thresholdIndex, window_median_length_index, onset  + round(window_length*2):offset - round(window_length*2)));
                fpRegion2 = squeeze(detection_res(loopIndex, JNRIndex, thresholdIndex, window_median_length_index, offset + round(window_length*2):end));
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

xAux = linspace(0, 1, numel(averageTPR(i, :, k)));
yAux = linspace(0, 1, numel(averageFPR(i, :, k)));
plot(xAux, yAux, '--');
box 'off'
grid on
xlabel('Probability of false alarm');
ylabel('Probability of detection');
legend('-20 dB', '-15 dB', '-10 dB', '-5 dB', '0 dB', 'Random guess', 'Location', 'southeast');


figure;
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
    hold on    
end

xAux = linspace(0, 1, numel(averageTPR(i, :, k)));
yAux = linspace(0, 1, numel(averageFPR(i, :, k)));
plot(xAux, yAux, '--');
box 'off'
grid on
xlabel('Probability of false alarm');
ylabel('Probability of detection');
legend('-20 dB', '-15 dB', '-10 dB', '-5 dB', '0 dB', 'Random guess', 'Location', 'southeast');


figure;
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
    hold on 
end

xAux = linspace(0, 1, numel(averageTPR(i, :, k)));
yAux = linspace(0, 1, numel(averageFPR(i, :, k)));
plot(xAux, yAux, '--');
box 'off'
grid on
xlabel('Probability of false alarm');
ylabel('Probability of detection');
legend('-20 dB', '-15 dB', '-10 dB', '-5 dB', '0 dB', 'Random guess', 'Location', 'southeast');


rmpath(['..' filesep '.' filesep 'data' filesep '07-24']);    
rmpath(['..' filesep '.' filesep 'Misc'])
