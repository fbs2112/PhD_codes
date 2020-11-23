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

load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_56.mat']);
load (['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'pfa_det_57.mat']);

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
dataPath = ['..' filesep '.' filesep 'figs' filesep '2020-11-05' filesep];

sizeTP = size(detection_res, 1);

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
p1 = plot(averageFpr.', averageTpr);
hold on;

thresholdTemplate = thresholdVector;
idxaux = find(~c);
if ~isempty(idxaux)
    idxProp = idxaux;
end

% load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_12.mat']);
% load (['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'pfa_det_13.mat']);
% tp = zeros(length(thresholdVector), sizeTP);
% fn = zeros(length(thresholdVector), sizeTP);
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


load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_48.mat']);
load (['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'pfa_det_49.mat']);
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
p2 = plot(averageFpr.', averageTpr);

thresholdFreq = thresholdVector;
idxaux = find(~c);
if ~isempty(idxaux)
    idxFreq = idxaux;
end


load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_52.mat']);
load (['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'pfa_det_53.mat']);

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
p3 = plot(averageFpr.', averageTpr);


thresholdTF = thresholdVector;
idxaux = find(~c);
if ~isempty(idxaux)
    idxTF = idxaux;
end

ylabel('Detection Probability');
xlabel('False Alarm Probability');
grid on
% legend('Template-based', 'Time-based', 'Frequency-based', 'TF-based', 'location', 'best');
% legend('Template-based', 'Frequency-based', 'TF-based', 'location', 'best');

colourPlot = get(gca,'colororder');

disp('--------------------------------------------------------------------')
% formatFig(gcf, [dataPath 'det_results_parkes_roc_new_train_test'], 'en', figProp);


%---------------------------Testing----------------------------------------

load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_58.mat']);
load (['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'pfa_det_59.mat']);

sizeTP2 = size(detection_res, 1);
tp = zeros(length(thresholdVector(idxProp)), sizeTP2);
fn = zeros(length(thresholdVector(idxProp)), sizeTP2);

for thresholdIndex = 1:length(thresholdVector(idxProp))
    x = squeeze(detection_res(:, idxProp(thresholdIndex)));
    tp(thresholdIndex, :) = x;
    fn(thresholdIndex, :) = size(x, 2) - tp(thresholdIndex, :);
end

tpr = squeeze(tp./(tp+fn));
averageTpr = mean(tpr, 2);

averageFpr = averageFpr(idxProp);

c = cfun(averageTpr, averageFpr.');

[cMinProp, ~] = min(c(:));
disp('-------------Testing-------------------------------------------')
disp('Proposed');
disp(['Pd: ' num2str(averageTpr)])
disp(['Pf: ' num2str(averageFpr)])
disp(['Cmin: ' num2str(cMinProp)]);

% load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_14.mat']);
% load (['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'pfa_det_15.mat']);
% 
% tp = zeros(length(thresholdVector(idxTime)), sizeTP2);
% fn = zeros(length(thresholdVector(idxTime)), sizeTP2);
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

load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_50.mat']);
load (['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'pfa_det_51.mat']);


tp = zeros(length(thresholdVector(idxFreq)), sizeTP2);
fn = zeros(length(thresholdVector(idxFreq)), sizeTP2);

for thresholdIndex = 1:length(thresholdVector(idxFreq))
    x = squeeze(detection_res(:, idxFreq(thresholdIndex)));
    tp(thresholdIndex, :) = x;
    fn(thresholdIndex, :) = size(x, 2) - tp(thresholdIndex, :);
end

tpr = squeeze(tp./(tp+fn));
averageTpr = mean(tpr, 2);
averageFpr = averageFpr(idxFreq);

c = cfun(averageTpr, averageFpr.');
[cMinFreq, ~] = min(c(:));

disp('Freq');
disp(['Pd: ' num2str(averageTpr)])
disp(['Pf: ' num2str(averageFpr)])
disp(['Cmin: ' num2str(cMinFreq)]);

load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_54.mat']);
load (['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'pfa_det_55.mat']);

tp = zeros(length(thresholdVector(idxTF)), sizeTP2);
fn = zeros(length(thresholdVector(idxTF)), sizeTP2);

for thresholdIndex = 1:length(thresholdVector(idxTF))
    x = squeeze(detection_res(:, idxTF(thresholdIndex)));
    tp(thresholdIndex, :) = x;
    fn(thresholdIndex, :) = size(x, 2) - tp(thresholdIndex, :);
end

tpr = squeeze(tp./(tp+fn));
averageTpr = mean(tpr, 2);
averageFpr = averageFpr(idxTF);

c = cfun(averageTpr, averageFpr.');
[cMinTF, ~] = min(c(:));

disp('TF');
disp(['Pd: ' num2str(averageTpr)])
disp(['Pf: ' num2str(averageFpr)])
disp(['Cmin: ' num2str(cMinTF)]);





load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_56.mat']);
load (['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'pfa_det_57.mat']);

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
dataPath = ['..' filesep '.' filesep 'figs' filesep '2020-11-05' filesep];

sizeTP = size(detection_res, 1);

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

axes('Position',[.53 .23 .34 .34])
box on
% figure;
semilogx(averageFpr.', averageTpr);
grid on

hold on;


load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_52.mat']);
load (['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'pfa_det_53.mat']);

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
% ylabel('Detection Probability');
% xlabel('False Alarm Probability');

semilogx(averageFpr.', averageTpr, '--', 'color', colourPlot(3,:));

xlim([1e-3 1e-1]);
legend([p1 p2 p3], 'Template-based', 'Frequency-based', 'TF-based', 'location', 'southwest');

% xlim([-12 -8]);
% ylim([0 0.3]);
% xticks([-12 -10 -8])
% xticklabels({'-12', '-10', '-8'});
% legend('Template-based', 'Time-based', 'Frequency-based', 'TF-based', 'location', 'best');
% legend('Template-based', 'TF-based', 'location', 'best');

% formatFig(gcf, [dataPath 'det_results_parkes_roc_new_train_test_log'], 'en', figProp);
formatFig(gcf, [dataPath 'det_results_parkes_roc_new_train_test'], 'en', figProp);
formatFig(gcf, [dataPath 'det_results_parkes_roc_new_train_test'], 'en', figProp);


rmpath(['..' filesep '.' filesep 'Misc'])