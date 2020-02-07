clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');

addpath(['..' filesep '.' filesep 'Misc'])

load(['.' filesep 'data' filesep 'detection01.mat']);

tpr = tpPR./(tpPR+fnPR);
averageTpr = squeeze(mean(tpr, 3));
stdTpr = std(tpr, [], 3);

fpr = tnPR./(fpPR + tnPR);

averageFpr = squeeze(mean(fpr, 3));
stdFpr = std(fpr, [], 3);

JNRVector = [-5 10 30 50];

for JNRIndex = 1:length(JNRVector)
    figure;

    for i = 1:3
%         plot(averageFpr(i,:), squeeze(averageTpr(i,JNRIndex,:)));
        loglog(averageFpr(i,:), squeeze(averageTpr(i,JNRIndex,:)));

        hold on;
    end
%     plot(linspace(0, 1, size(averageFpr, 2)), linspace(0, 1, size(averageFpr, 2)));
%     loglog(logspace(0, 1, size(averageFpr, 2)), logspace(0, 1, size(averageFpr, 2)));

    ylabel('Probability of detection');
    xlabel('Probability of false alarm');
    grid on;
    title(['JNR = ' num2str(JNRVector(JNRIndex))]);
    legend('NMF-Based', 'Pulse Blanking', 'No mitigation', 'Random guess');
end

%--------------------------------------------------------------------------

tpr = tpPRB./(tpPRB+fnPRB);
averageTpr = squeeze(mean(tpr, 3));
stdTpr = std(tpr, [], 3);

fpr = tnPRB./(fpPRB + tnPRB);

averageFpr = squeeze(mean(fpr, 3));
stdFpr = std(fpr, [], 3);

JNRVector = [-5 10 30 50];

for JNRIndex = 1:length(JNRVector)
    figure;

    for i = 1:3
        plot(averageFpr(i,:), squeeze(averageTpr(i,JNRIndex,:)));
        hold on;
        plot(linspace(0, 1, numel(averageFpr)), linspace(0, 1, numel(averageFpr)));
        ylabel('Probability of detection');
        xlabel('Probability of false alarm');
        grid on;
        title(['JNR = ' num2str(JNRVector(JNRIndex))]);
    end
end
%--------------------------------------------------------------------------

tpr = tpPAPR./(tpPAPR+fnPAPR);
averageTpr = squeeze(mean(tpr, 3));
stdTpr = std(tpr, [], 3);

fpr = tnPAPR./(fpPAPR + tnPAPR);

averageFpr = squeeze(mean(fpr, 3));
stdFpr = std(fpr, [], 3);

JNRVector = [-5 10 30 50];

for JNRIndex = 1:length(JNRVector)
    for i = 1:3
        figure;
        plot(averageFpr(i,:), squeeze(averageTpr(i,JNRIndex,:)));
        hold on;
        plot(linspace(0, 1, numel(averageFpr)), linspace(0, 1, numel(averageFpr)));
        ylabel('Probability of detection');
        xlabel('Probability of false alarm');
        grid on;
        title(['JNR = ' num2str(JNRVector(JNRIndex))]);
    end
end

rmpath(['..' filesep '.' filesep 'Misc']);