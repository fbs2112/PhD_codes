clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '.' filesep 'Misc'])

frequency = '1152';
fileName = '2018-09-01-09_52_55_0000020535312384.000000';

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
dataPath = ['..' filesep '.' filesep 'figs' filesep '2020-11-05' filesep];

load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_58.mat']);
figure;
% histogram(feature, 10);
[f,xi] = ksdensity(feature);
plot(xi, f);
hold on;

load (['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_59.mat']);
[f,xi] = ksdensity(feature);
plot(xi, f);
% histogram(feature, 20);
xline(0.3595, '--r', 'linewidth', 1.5)
legend('ADS-B', 'No ADS-B', '$\bar{\gamma}$');
xlabel('$s_{\mathrm{max}}$')
ylabel('Occurrences');
% formatFig(gcf, [dataPath 'hist_template'], 'en', figProp);

clearvars feature featurePfa;
load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_54.mat']);
figure;
% [f,xi] = ksdensity(feature(feature < 2000));
[f,xi] = ksdensity(feature);
plot(xi, f);
% histogram(feature(feature < 2000), 10);
hold on;

load (['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_55.mat']);
[f,xi] = ksdensity(featurePfa);
plot(xi, f);
% histogram(featurePfa, 20);
xline(764.5, '--r', 'linewidth', 1.5)
legend('ADS-B', 'No ADS-B', '$\bar{\gamma}$');
% xlim([500 2000]);
xlabel('$X_{\mathrm{TF_{max}}}$')
ylabel('Occurrences');
% formatFig(gcf, [dataPath 'hist_tf'], 'en', figProp);

clearvars feature featurePfa;
load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_50.mat']);
figure;
[f,xi] = ksdensity(feature);
plot(xi, f);
% histogram(feature, 10);
hold on;

load (['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_51.mat']);
[f,xi] = ksdensity(featurePfa);
plot(xi, f);
% histogram(featurePfa, 20);
xline(17372, '--r', 'linewidth', 1.5)
legend('ADS-B', 'No ADS-B', '$\bar{\gamma}$');
xlabel('$X_{\mathrm{f_{max}}}$')
ylabel('Occurrences');
% formatFig(gcf, [dataPath 'hist_freq'], 'en', figProp);


rmpath(['..' filesep '.' filesep 'Misc'])