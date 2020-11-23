clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '.' filesep 'Misc'])

cfun = @(tpr, fpr) sqrt(fpr.^2 + (1-tpr).^2);

frequency = '1152';
fileName = '2018-09-01-09_52_55_0000020535312384.000000';

load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_46.mat']);
load (['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'pfa_det_47.mat']);

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
dataPath = ['..' filesep '.' filesep 'figs' filesep '2020-11-05' filesep];

% thresholdVector = 0.4315;
sizeTP = 105;

tp = zeros(length(thresholdVector), sizeTP);
fn = zeros(length(thresholdVector), sizeTP);

for thresholdIndex = 1:length(thresholdVector)
    x = squeeze(detection_res(:, thresholdIndex));
    tp(thresholdIndex, :) = x;
    fn(thresholdIndex, :) = size(x, 2) - tp(thresholdIndex, :);
end

tpr = squeeze(tp./(tp+fn));
averageTpr = mean(tpr, 2);
stdTpr = std(tpr, [], 2);

c = cfun(averageTpr, averageFpr.');
[cMinProp, idxProp] = min(c(:));
disp('-------------Training-------------------------------------------')
disp('Proposed');
disp(['Pd: ' num2str(averageTpr(idxProp))])
disp(['Pf: ' num2str(averageFpr(idxProp))])
disp(['Cmin: ' num2str(cMinProp)]);

figure;
plot(averageFpr.', averageTpr);
hold on;

thresholdTemplate = thresholdVector;
idxaux = find(~c);
if ~isempty(idxaux)
    idxProp = idxaux;
end

% load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_14.mat']);
% load (['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'pfa_det_15.mat']);
% tp = zeros(length(thresholdVector), 32);
% fn = zeros(length(thresholdVector), 32);
% for thresholdIndex = 1:length(thresholdVector)
% %     thresholdIndex
%     x = squeeze(detection_res(:, thresholdIndex));
%     tp(thresholdIndex, :) = x;
%     fn(thresholdIndex, :) = 1 - tp(thresholdIndex, :);
% end
% 
% tpr = squeeze(tp./(tp+fn));
% averageTpr = mean(tpr, 2);
% c = cfun(averageTpr, averageFpr.');
% [cMinTime, idxTime] = min(c(:));
% 
% disp('Time');
% disp(['Pd: ' num2str(averageTpr(idxTime))])
% disp(['Pf: ' num2str(averageFpr(idxTime))])
% disp(['Cmin: ' num2str(cMinTime)]);
% 
% plot(averageFpr.', averageTpr);
% hold on
% 
% thresholdTime = thresholdVector;
% idxaux = find(~c);
% if ~isempty(idxaux)
%     idxTime = idxaux;
% end


load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_38.mat']);
load (['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'pfa_det_39.mat']);
tp = zeros(length(thresholdVector), sizeTP);
fn = zeros(length(thresholdVector), sizeTP);
for thresholdIndex = 1:length(thresholdVector)
    x = squeeze(detection_res(:, thresholdIndex));
    tp(thresholdIndex, :) = x;
    fn(thresholdIndex, :) = size(x, 2) - tp(thresholdIndex, :);
end

tpr = squeeze(tp./(tp+fn));
averageTpr = mean(tpr, 2);
c = cfun(averageTpr, averageFpr.');
[cMinFreq, idxFreq] = min(c(:));
disp('Freq');
disp(['Pd: ' num2str(averageTpr(idxFreq))])
disp(['Pf: ' num2str(averageFpr(idxFreq))])
disp(['Cmin: ' num2str(cMinFreq)]);
plot(averageFpr.', averageTpr);

thresholdFreq = thresholdVector;
idxaux = find(~c);
if ~isempty(idxaux)
    idxFreq = idxaux;
end


load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_42.mat']);
load (['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'pfa_det_43.mat']);

tp = zeros(length(thresholdVector), sizeTP);
fn = zeros(length(thresholdVector), sizeTP);
for thresholdIndex = 1:length(thresholdVector)
    x = squeeze(detection_res(:, thresholdIndex));
    tp(thresholdIndex, :) = x;
    fn(thresholdIndex, :) = size(x, 2) - tp(thresholdIndex, :);
end

tpr = squeeze(tp./(tp+fn));
averageTpr = mean(tpr, 2);
c = cfun(averageTpr, averageFpr.');
[cMinTF, idxTF] = min(c(:));
disp('TF');
disp(['Pd: ' num2str(averageTpr(idxTF))])
disp(['Pf: ' num2str(averageFpr(idxTF))])
disp(['Cmin: ' num2str(cMinTF)]);
plot(averageFpr.', averageTpr);


thresholdTF = thresholdVector;
idxaux = find(~c);
if ~isempty(idxaux)
    idxTF = idxaux;
end

ylabel('True Positive Rate');
xlabel('False Positive Rate');
grid on
legend('Template-based', 'Frequency-based', 'TF-based', 'location', 'best');
% legend('Template-based', 'Time-based', 'Frequency-based', 'TF-based', 'location', 'best');


disp('--------------------------------------------------------------------')
% formatFig(gcf, [dataPath 'det_results_parkes_roc_new'], 'en', figProp);


%---------------------------Testing----------------------------------------

load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_44.mat']);
load (['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'pfa_det_45.mat']);

sizeTP = 19;

tp = zeros(length(thresholdVector(idxProp)), sizeTP);
fn = zeros(length(thresholdVector(idxProp)), sizeTP);

for thresholdIndex = 1:length(thresholdVector(idxProp))
    x = squeeze(detection_res(:, idxProp(thresholdIndex)));
    tp(thresholdIndex, :) = x;
    fn(thresholdIndex, :) = size(x, 2) - tp(thresholdIndex, :);
end

tpr = squeeze(tp./(tp+fn));
averageTpr = mean(tpr, 2);

averageFpr = averageFpr(idxProp);

c = cfun(averageTpr, averageFpr.');

[cMinProp, idxProp] = min(c(:));
disp('-------------Testing-------------------------------------------')
disp('Proposed');
disp(['Pd: ' num2str(averageTpr(idxProp))])
disp(['Pf: ' num2str(averageFpr(idxProp))])
disp(['Cmin: ' num2str(cMinProp)]);

% load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_12.mat']);
% load (['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'pfa_det_13.mat']);
% 
% tp = zeros(length(thresholdVector(idxTime)), 10);
% fn = zeros(length(thresholdVector(idxTime)), 10);
% 
% for thresholdIndex = 1:length(thresholdVector(idxTime))
%     x = squeeze(detection_res(:, idxTime(thresholdIndex)));
%     tp(thresholdIndex, :) = x;
%     fn(thresholdIndex, :) = size(x, 2) - tp(thresholdIndex, :);
% end
% 
% tpr = squeeze(tp./(tp+fn));
% averageTpr = mean(tpr, 2);
% averageFpr = averageFpr(idxTime);
% 
% c = cfun(averageTpr, averageFpr.');
% [cMinTime, idxTime] = min(c(:));
% disp('Time');
% disp(['Pd: ' num2str(averageTpr(idxTime))])
% disp(['Pf: ' num2str(averageFpr(idxTime))])
% disp(['Cmin: ' num2str(cMinTime)]);

load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_36.mat']);
load (['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'pfa_det_37.mat']);


tp = zeros(length(thresholdVector(idxFreq)), sizeTP);
fn = zeros(length(thresholdVector(idxFreq)), sizeTP);

for thresholdIndex = 1:length(thresholdVector(idxFreq))
    x = squeeze(detection_res(:, idxFreq(thresholdIndex)));
    tp(thresholdIndex, :) = x;
    fn(thresholdIndex, :) = size(x, 2) - tp(thresholdIndex, :);
end

tpr = squeeze(tp./(tp+fn));
averageTpr = mean(tpr, 2);
averageFpr = averageFpr(idxFreq);

c = cfun(averageTpr, averageFpr.');
[cMinFreq, idxFreq] = min(c(:));


disp('Freq');
disp(['Pd: ' num2str(averageTpr(idxFreq))])
disp(['Pf: ' num2str(averageFpr(idxFreq))])
disp(['Cmin: ' num2str(cMinFreq)]);

load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_40.mat']);
load (['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'pfa_det_41.mat']);

tp = zeros(length(thresholdVector(idxTF)), sizeTP);
fn = zeros(length(thresholdVector(idxTF)), sizeTP);

for thresholdIndex = 1:length(thresholdVector(idxTF))
    x = squeeze(detection_res(:, idxTF(thresholdIndex)));
    tp(thresholdIndex, :) = x;
    fn(thresholdIndex, :) = size(x, 2) - tp(thresholdIndex, :);
end

tpr = squeeze(tp./(tp+fn));
averageTpr = mean(tpr, 2);
averageFpr = averageFpr(idxTF);

c = cfun(averageTpr, averageFpr.');
[cMinTF, idxTF] = min(c(:));

disp('TF');
disp(['Pd: ' num2str(averageTpr(idxTF))])
disp(['Pf: ' num2str(averageFpr(idxTF))])
disp(['Cmin: ' num2str(cMinTF)]);

rmpath(['..' filesep '.' filesep 'Misc'])