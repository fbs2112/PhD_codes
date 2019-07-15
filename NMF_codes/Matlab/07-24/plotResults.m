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