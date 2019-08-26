%Probability of false alarm evaluation for pai code

clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '..' filesep '.' filesep 'Misc'])
addpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'pai_results']);  

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);

load results13.mat;

WinLBlock = [3 19];
PfaVector = logspace(-5, 0, 17);
monteCarloLoops = 1000;

for k = 1:length(WinLBlock)
    x = squeeze(detection_res_cell{k, 1});
    idx = randperm(size(x, 2), 1);
    fp = zeros(length(PfaVector), monteCarloLoops);
    tn = zeros(length(PfaVector), monteCarloLoops);
    for j = 1:length(PfaVector)
        for i = 1:monteCarloLoops
%             idx = randperm(size(x, 2), 1);
            fp(j,i) = sum(x(i,idx,j));
            tn(j,i) = 1 - fp(j,i);
        end
    end
    fpr(k,:,:) = fp./(fp+tn);
end

averageFpr = squeeze(mean(fpr, 3));
stdFpr = squeeze(std(fpr, [], 3));

figure;
for i = 1:size(averageFpr, 1)
    plot(PfaVector, averageFpr(i,:))
    hold on;
end

ylabel('Probability of false alarm');
xlabel('$\bar{\gamma}$');
xlim([min(PfaVector) max(PfaVector)]);
legend('$L_{\mathrm{STFT}} = 3$', '$L_{\mathrm{STFT}} = 19$');
grid on;

figure;
for i = 1:size(averageFpr, 1)
    loglog(PfaVector, averageFpr(i,:))
    hold on;
end
plot(logspace(-5, 0, numel(PfaVector)), logspace(-5, 0, numel(PfaVector)));

ylabel('Probability of false alarm');
xlabel('$\bar{\gamma}$');
xlim([min(PfaVector) max(PfaVector)]);
ylim([1e-5 1e-0]);
legend('$L_{\mathrm{STFT}} = 3$', '$L_{\mathrm{STFT}} = 19$', 'Random guess');
grid on;

save(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'pai_results' filesep 'pfa_data_2.mat'], 'averageFpr', 'stdFpr')

rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])
rmpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'pai_results']);  