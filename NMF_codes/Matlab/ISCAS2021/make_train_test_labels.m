clear;
clc;
close all;

frequency = '1152';
% fileName = '2018-09-01-09_52_55_0000020535312384.000000'; %train
fileName = '2018-09-01-09_52_55_0000000000000000.000000'; %test

load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep fileName '_labels.mat']);
save(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep fileName '_labels_bkp.mat'], 'interferenceDetFlag');

falseLabels = find(interferenceDetFlag == 0 | interferenceDetFlag == 1 | interferenceDetFlag == 4);
trueLabels = find(interferenceDetFlag == 2 | interferenceDetFlag == 3);

falseLabelsLength = length(falseLabels);
trueLabelsLength = length(trueLabels);

proportionTrain = 0.6;
proportionTest = 0.4;

falseLabelsLengthTrain = round(falseLabelsLength * proportionTrain);
trueLabelsLengthTrain = round(trueLabelsLength * proportionTrain);

falseLabelsLengthTest = falseLabelsLength - falseLabelsLengthTrain;
trueLabelsLengthTest = trueLabelsLength - trueLabelsLengthTrain;

%Label creation

aux = randperm(falseLabelsLength);

falseLabelsTrain = falseLabels(aux(1:falseLabelsLengthTrain));
falseLabelsTest = falseLabels(aux(falseLabelsLengthTrain+1:end));

aux = randperm(trueLabelsLength);

trueLabelsTrain = trueLabels(aux(1:trueLabelsLengthTrain));
trueLabelsTest = trueLabels(aux(trueLabelsLengthTrain+1:end));

intersect(falseLabelsTrain, falseLabelsTest)
intersect(falseLabelsTrain, trueLabelsTrain)
intersect(falseLabelsTrain, trueLabelsTest)
intersect(falseLabelsTest, trueLabelsTrain)
intersect(falseLabelsTest, trueLabelsTest)
intersect(trueLabelsTrain, trueLabelsTest)

if isempty(intersect(falseLabelsTrain, falseLabelsTest)) && isempty(intersect(trueLabelsTrain, trueLabelsTest))
   save(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep fileName '_labels_train_test.mat'],...
       'falseLabelsTrain', 'falseLabelsTest', 'trueLabelsTrain', 'trueLabelsTest');
end