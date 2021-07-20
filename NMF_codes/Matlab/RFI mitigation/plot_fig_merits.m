clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');

addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'Misc'])
addpath(['..' filesep 'Sigtools' filesep])

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
dataPath = ['..' filesep 'figs' filesep '2020-06-10' filesep];
figProp = struct('size' , fontsize , 'font' , fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);

load(['.' filesep 'data' filesep 'fig_merit_1.mat']);

bandwidthVector = [2 8 14];
SNR = -20;
JNRVector = (0:5:30) - SNR;

averageGenSNR = mean(generalisedSNR, 4);

M = 20;
for bandwidthIndex = 1:length(bandwidthVector)
    for JNRIndex = 1:length(JNRVector)
        for i = 1:7
            corrBuffer = buffer(corrOut(i,JNRIndex,bandwidthIndex,:), M, 0, 'nodelay');
            Pn = sum(real(corrBuffer), 1).^2 + sum(imag(corrBuffer), 1).^2;
            Pw = sum(real(corrBuffer).^2 + imag(corrBuffer).^2, 1);
            s(i,JNRIndex,bandwidthIndex) = mean(Pn ./ Pw);
            stdS(i,JNRIndex,bandwidthIndex) = std(Pn ./ Pw);
            if s(i,JNRIndex,bandwidthIndex) < 1
                s(i,JNRIndex,bandwidthIndex) = 1;
            end
            C_N0(i,JNRIndex,bandwidthIndex) = 10*log10((1/1e-3) * ((s(i,JNRIndex,bandwidthIndex) - 1))/((M - s(i,JNRIndex,bandwidthIndex))));
        end
    end
end


for bandwidthIndex = 1:length(bandwidthVector)
    figure;
    plot(JNRVector, pow2db(averageGenSNR(1,:,bandwidthIndex)), '-o', 'Color', [0 0.4470 0.7410]);
    hold on;
    plot(JNRVector, pow2db(averageGenSNR(2,:,bandwidthIndex)), '--+', 'Color', [0.850  0.3250 0.0980]);
    plot(JNRVector, pow2db(averageGenSNR(3,:,bandwidthIndex)), '--*', 'Color', [0.9290 0.6940 0.1250]);
    plot(JNRVector, pow2db(averageGenSNR(4,:,bandwidthIndex)), '-^', 'Color', [0.4940 0.1840 0.5560]);
    plot(JNRVector, pow2db(averageGenSNR(5,:,bandwidthIndex)), '-x', 'Color', [0.4660 0.6740 0.1880]);
    plot(JNRVector, pow2db(averageGenSNR(6,:,bandwidthIndex)), '-s', 'Color', [0.3010 0.7450 0.9330]);
    plot(JNRVector, pow2db(averageGenSNR(7,:,bandwidthIndex)), '-d', 'Color', [0.6350 0.0780 0.1840]);
    grid on;
    ylim([-10 50]);
    ylabel('Generalised SNR [dB]')
    xlabel('JSR [dB]');
    legend('NMF-stft', 'Kalman-Perf', 'Kalman-Est.', 'NMF-fsst (Kaiser)', 'NMF-stft (Kaiser)','NMFSB-fsst (Kaiser)', 'NMFSB-stft (Kaiser)', 'location', 'best');
%     formatFig(gcf, [dataPath 'av_snr_' num2str(bandwidthIndex)], 'en', figProp);
end

for bandwidthIndex = 1:length(bandwidthVector)
    figure;
    plot(JNRVector, squeeze(C_N0(1,:,bandwidthIndex)), '-o', 'Color', [0 0.4470 0.7410]);
    hold on;
    plot(JNRVector, squeeze(C_N0(2,:,bandwidthIndex)), '--+', 'Color', [0.8500 0.3250 0.0980]);
    plot(JNRVector, squeeze(C_N0(3,:,bandwidthIndex)), '--*', 'Color', [0.9290 0.6940 0.1250]);
    plot(JNRVector, squeeze(C_N0(4,:,bandwidthIndex)), '-^', 'Color', [0.4940 0.1840 0.5560]);
    plot(JNRVector, squeeze(C_N0(5,:,bandwidthIndex)), '-x', 'Color', [0.4660 0.6740 0.1880]);
    plot(JNRVector, squeeze(C_N0(6,:,bandwidthIndex)), '-s', 'Color', [0.3010 0.7450 0.9330]);
    plot(JNRVector, squeeze(C_N0(7,:,bandwidthIndex)), '-d', 'Color', [0.6350 0.0780 0.1840]);
    grid on;
    ylabel('$C/N_0$ [dB]')
    xlabel('JSR [dB]');
    legend('NMF-stft', 'Kalman-Perf', 'Kalman-Est.', 'NMF-fsst (Kaiser)', 'NMF-stft (Kaiser)','NMFSB-fsst (Kaiser)', 'NMFSB-stft (Kaiser)', 'location', 'best');
end


rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'Misc'])
rmpath(['..' filesep 'signalsGeneration' filesep]);