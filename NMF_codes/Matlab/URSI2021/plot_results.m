clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');

addpath(['..' filesep 'Misc'])

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
dataPath = ['..' filesep 'figs' filesep 'thesis' filesep];
figProp = struct('size' , fontsize , 'font' , fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);

load(['.' filesep 'data' filesep 'results_10_11.mat']);
load(['.' filesep 'data' filesep 'discarded_10_11.mat']);

JNRVector = 30:5:60;
markerMatrix = {'o','+','*','^','s','v','p'};

figure;
plot(JNRVector, averagePPR(1,:), 'Color', [0 0.4470 0.7410]);
hold on;
plot(JNRVector, averagePPR(2,:), 'Color', [0.8500 0.3250 0.0980]);
plot(JNRVector, averagePPR(3,:), 'Color', [0.9290 0.6940 0.1250]);
plot(JNRVector, averagePPR(4,:), 'Color', [0.4940, 0.1840, 0.5560]);
plot(JNRVector, averagePPR(5,:), 'Color', [0.4660 0.6740 0.1880]);
grid on;
ylabel('Peak-to-peak Ratio')
xlabel('JSR [dB]');
legend('NMF-256', 'NMF-1024', 'PB-time', 'PB-freq', 'No mitigation', 'location', 'best');

figure;
plot(JNRVector, 10*log10(averageGenSNR(1,:)), 'Marker', markerMatrix{1}, 'Color', [0 0.4470 0.7410], 'Markersize', 12);
hold on;
plot(JNRVector, 10*log10(averageGenSNR(2,:)), 'Marker', markerMatrix{2}, 'Color', [0.8500 0.3250 0.0980], 'Markersize', 12);
plot(JNRVector, 10*log10(averageGenSNR(3,:)), 'Marker', markerMatrix{3}, 'Color', [0.9290 0.6940 0.1250], 'Markersize', 12);
plot(JNRVector, 10*log10(averageGenSNR(4,:)), 'Marker', markerMatrix{4}, 'Color', [0.4940 0.1840 0.5560], 'Markersize', 12);
plot(JNRVector, 10*log10(averageGenSNR(5,:)), 'Marker', markerMatrix{5}, 'Color', [0.4660 0.6740 0.1880], 'Markersize', 12);
grid on;
ylabel('GSNR [dB]')
xlabel('JSR [dB]');
legend('NMF-256', 'NMF-1024', 'PB-time', 'PB-freq', 'No mitigation', 'location', 'best');
formatFig(gcf, [dataPath  'gen_snr'], 'en', figProp);


figure;
plot(JNRVector, 10*log10(averageGenSNRMaxPeak(1,:)), 'Color', [0 0.4470 0.7410]);
hold on;
plot(JNRVector, 10*log10(averageGenSNRMaxPeak(2,:)), 'Color', [0.8500 0.3250 0.0980]);
plot(JNRVector, 10*log10(averageGenSNRMaxPeak(3,:)), 'Color', [0.9290 0.6940 0.1250]);
plot(JNRVector, 10*log10(averageGenSNRMaxPeak(4,:)), 'Color', [0.4940, 0.1840, 0.5560]);
plot(JNRVector, 10*log10(averageGenSNRMaxPeak(5,:)), 'Color', [0.4660 0.6740 0.1880]);

grid on;
ylabel('Generalized SNR Max Peak [dB]')
xlabel('JSR [dB]');
legend('NMF-256', 'NMF-1024', 'PB-time', 'PB-freq', 'No mitigation', 'location', 'best');

figure;
plot(JNRVector, peakToPeakRatio(1,:), 'Color', [0 0.4470 0.7410]);
hold on;
plot(JNRVector, peakToPeakRatio(2,:), 'Color', [0.8500 0.3250 0.0980]);
plot(JNRVector, peakToPeakRatio(3,:), 'Color', [0.9290 0.6940 0.1250]);
plot(JNRVector, peakToPeakRatio(4,:), 'Color', [0.4940, 0.1840, 0.5560]);
plot(JNRVector, peakToPeakRatio(5,:), 'Color', [0.4660 0.6740 0.1880]);

ylabel('Peak-to-peak Ratio')        
xlabel('JSR [dB]');
legend('NMF-256', 'NMF-1024', 'PB-Time', 'PB-Freq', 'No mitigation', 'location', 'best');

figure;
bar(JNRVector, C_N0.');
ylabel('C/N$_0$ [dB]')
xlabel('JSR [dB]');
legend('NMF-256', 'NMF-1024', 'PB-Time', 'PB-Freq', 'No mitigation', 'location', 'best');
grid on;
formatFig(gcf, [dataPath  'c_no'], 'en', figProp);
formatFig(gcf, [dataPath  'c_no'], 'en', figProp);

meanNumberOfZeroedSamples = mean(numberOfZeroedSamples, 2);
meanNumberOfZeroedSamplesFreq = mean(numberOfZeroedSamplesFreq, 2);

aux = [meanNumberOfZeroedSamples meanNumberOfZeroedSamplesFreq];

figure;
b = bar(JNRVector, (aux ./ 50000)*100);
b(1).FaceColor = [0.9290 0.6940 0.1250];
b(2).FaceColor = [0.4940, 0.1840, 0.5560];

ylabel('Discarded samples [\%]')
xlabel('JSR [dB]');
legend('PB-time', 'PB-freq', 'location', 'best');
grid on;
formatFig(gcf, [dataPath  'discard_samples'], 'en', figProp);
formatFig(gcf, [dataPath  'discard_samples'], 'en', figProp);

rmpath(['..' filesep 'Misc'])