%%
%Plot results for CW considering the new approach of using a large median
%window

clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '..' filesep '.' filesep 'Misc'])
<<<<<<< Updated upstream
addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
    'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
    filesep 'pfa_results' filesep]);
addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
    'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
    filesep]);
=======
if isunix
    addpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/my_results/pfa_results/');
    addpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/my_results/');
    
else
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep 'pfa_results' filesep])
>>>>>>> Stashed changes
    

load results_det_19.mat;
load resultspfa14.mat;

averageFprDot = averageFpr;

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
dataPath = ['..' filesep '..' filesep '.' filesep 'figs' filesep '08-12' filesep];

monteCarloLoops = 100;
thresholdVector = 0:0.005:0.2;
periodVector = 0;
bandwidthVector = 0;
JNRVector = -25:0;

for JNRIndex = 1:length(JNRVector)
    for bandwidthIndex = 1:length(bandwidthVector)
        for periodIndex = 1:length(periodVector)
            for thresholdIndex = 1:length(thresholdVector)
                x = squeeze(detection_res(:, bandwidthIndex, periodIndex, JNRIndex, thresholdIndex));
                tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :) = x;
                fn(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :) = ...
                    size(x, 2) - tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :);
            end
        end
    end
end

tpr = squeeze(tp./(tp+fn));

averageTpr = mean(tpr, 3);
stdTpr = std(tpr, [], 3);

cfun = @(tpr, fpr) sqrt(fpr.^2 + (1-tpr).^2);

for i = 1:length(JNRVector)
    c = cfun(averageTpr(i,:), averageFpr);
    cMin(i) = min(c(:));
end

<<<<<<< Updated upstream

rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
    'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
    filesep 'pfa_results' filesep]);
rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
    'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
    filesep]);

=======
if isunix
    rmpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/my_results/pfa_results/');
    rmpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/my_results/');
    
else
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep 'pfa_results' filesep])
    
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep]);
end  
>>>>>>> Stashed changes


%-------------------------------Pai's results------------------------------

<<<<<<< Updated upstream
addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
    'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
    filesep 'pfa_results' filesep]);

addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
    'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
    filesep]);
=======
if isunix
    addpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/pai_results/pfa_results/');
    addpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/pai_results/');
    
else
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep 'pfa_results' filesep])
>>>>>>> Stashed changes
    
load resultspfa3.mat;
load results_det_7.mat;

averageFprPai = averageFpr;
JNRVector = -25:0;
PfaVector = logspace(-12, -2, 41);

for JNRIndex = 1:length(JNRVector)
    
    for i = 1:length(PfaVector)
        tpRegion = squeeze(detection_res(JNRIndex,:,:,i));
        tpPai(JNRIndex,i,:) = any(tpRegion, 2);
        fnPai(JNRIndex,i,:) = 1 - tpPai(JNRIndex,i,:);
    end

end

tprPai = tpPai./(tpPai+fnPai);
averageTprPai = mean(tprPai, 3);
stdTprPai = std(tprPai, [], 3);

for i = 1:length(JNRVector)
    c = cfun(averageTprPai(i,:), averageFpr(1,:));
    [cMinPai1(i), idx1(i)] = min(c(:));
end

load results_det_6.mat;

for JNRIndex = 1:length(JNRVector)
    
    for i = 1:length(PfaVector)
        tpRegion = squeeze(detection_res(JNRIndex,:,:,i));
        tpPai2(JNRIndex,i,:) = any(tpRegion, 2);
        fnPai2(JNRIndex,i,:) = 1 - tpPai2(JNRIndex,i,:);
    end

end

tprPai2 = tpPai2./(tpPai2+fnPai2);
averageTprPai2 = mean(tprPai2, 3);
stdTprPai2 = std(tprPai2, [], 3);

for i = 1:length(JNRVector)
    c = cfun(averageTprPai2(i,:), averageFpr(2,:));
    [cMinPai2(i), idx2(i)] = min(c(:));
end

<<<<<<< Updated upstream
load resultspfa4.mat;
load results_det_11.mat;
=======

load results_det_11.mat;
load resultspfa4.mat;
averageFprPai3 = averageFpr;
>>>>>>> Stashed changes

for JNRIndex = 1:length(JNRVector)
    
    for i = 1:length(PfaVector)
        tpRegion = squeeze(detection_res(JNRIndex,:,:,i));
        tpPai3(JNRIndex,i,:) = any(tpRegion, 2);
        fnPai3(JNRIndex,i,:) = 1 - tpPai3(JNRIndex,i,:);
    end

end

tprPai3 = tpPai3./(tpPai3+fnPai3);
averageTprPai3 = mean(tprPai3, 3);
stdTprPai3 = std(tprPai3, [], 3);

for i = 1:length(JNRVector)
<<<<<<< Updated upstream
    c = cfun(averageTprPai3(i,:), averageFpr);
    [cMinPai3(i), idx3(i)] = min(c(:));
end

=======
    c = cfun(averageTprPai3(i,:), averageFprPai3);
    [cMinPai3(i), idx3(i)] = min(c(:));
end



>>>>>>> Stashed changes
figure;
plot(JNRVector, cMin)
hold on;
plot(JNRVector, cMinPai1)
plot(JNRVector, cMinPai3)
plot(JNRVector, cMinPai2)
ylim([0 1/sqrt(2)])
grid on;
ylabel('C$_{\mathrm{min}}$');
xlabel('JNR [dB]');
legend('Dot', 'Statistical, $L = 3$', 'Statistical, $L = 11$', 'Statistical, $L = 19$');

formatFig(gcf, [dataPath  'cmin_cw'], 'en', figProp);

figure;
plot(JNRVector, cMin)
hold on;
plot(JNRVector, cMinPai1)
plot(JNRVector, cMinPai3)
plot(JNRVector, cMinPai2)
ylim([0 1/sqrt(2)])
xlim([-17 -12]);
grid on;
ylabel('C$_{\mathrm{min}}$');
xlabel('JNR [dB]');
legend('Dot', 'Statistical, $L = 3$', 'Statistical, $L = 11$', 'Statistical, $L = 19$');

formatFig(gcf, [dataPath  'cmin_cw_zoom'], 'en', figProp);

 for i = 1:length(JNRVector)
    figure;
    plot(averageFprDot, averageTpr(i,:));
    hold on;
    plot(averageFprPai(1,:), averageTprPai(i,:));
    plot(averageFprPai3, averageTprPai3(i,:));
    plot(averageFprPai(2,:), averageTprPai2(i,:));
    plot(linspace(0, 1, numel(averageFpr)), linspace(0, 1, numel(averageFpr)));
    legend('Dot', 'Statistical, $L = 3$', 'Statistical, $L = 11$', 'Statistical, $L = 19$', 'Random guess', 'location', 'southeast');
    ylabel('Probability of detection');
    xlabel('Probability of false alarm');
    grid on;
    formatFig(gcf, [dataPath  'roc_cw_all_' num2str(JNRVector(i))], 'en', figProp);
end

rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])
<<<<<<< Updated upstream
rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
    'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
    filesep 'pfa_results' filesep]);
=======
if isunix
    rmpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/pai_results/pfa_results/');
    rmpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/pai_results/');
    
else
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep 'pfa_results' filesep])
    
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep]);
end  
>>>>>>> Stashed changes

rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
    'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
    filesep]);
    
%%
%Plot results for chirp for different bandwidths and periods (my technique)

clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '..' filesep '.' filesep 'Misc'])
<<<<<<< Updated upstream
addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
    'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
    filesep 'pfa_results' filesep]);
addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
    'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
    filesep]);
=======
if isunix
    addpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/my_results/pfa_results/');
    addpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/my_results/');
    
else
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep 'pfa_results' filesep])
    
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep]);
end  

>>>>>>> Stashed changes

load results_det_23.mat;
load resultspfa14.mat;

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
dataPath = ['..' filesep '..' filesep '.' filesep 'figs' filesep '08-12' filesep];

thresholdVector = 0:0.005:0.2;
bandwidthVector = (2e6:3e6:14e6)/1e6;
periodVector = ((8.62e-6:1.48e-6:18.97e-6)*1e6);
JNRVector = -25:0;

for JNRIndex = 1:length(JNRVector)
    for bandwidthIndex = 1:length(bandwidthVector)
        for periodIndex = 1:length(periodVector)
            for thresholdIndex = 1:length(thresholdVector)
                x = squeeze(detection_res(:, bandwidthIndex, periodIndex, JNRIndex, thresholdIndex));
                tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :) = x;
                fn(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :) = ...
                    size(x, 2) - tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :);
            end
        end
    end
end

tpr = squeeze(tp./(tp+fn));
averageTpr = mean(tpr, 5);
stdTpr = std(tpr, [], 5);

cfun = @(tpr, fpr) sqrt(fpr.^2 + (1-tpr).^2);

for k = 1:length(bandwidthVector)
    for j = 1:length(periodVector)
        for i = 1:length(JNRVector)
            c = cfun(squeeze(averageTpr(k,j,i,:)).', averageFpr);
            [cMin(k,j,i), idx] = min(c(:));
        end
    end
end

for k = 1:length(bandwidthVector)
    figure;
    
    for j = 1:length(periodVector)
        plot(JNRVector, squeeze(cMin(k,j,:)));
        hold on;
    end
    ylabel('C$_{\mathrm{min}}$');
    xlabel('JNR [dB]');
    legend(['$T$ = ' num2str(periodVector(1), '%.2f') ' $\mu$s'], ['$T$ = ' num2str(periodVector(2), '%.2f') ' $\mu$s'], ['$T$ = ' num2str(periodVector(3), '%.2f') ' $\mu$s'],...
        ['$T$ = ' num2str(periodVector(4), '%.2f') ' $\mu$s'], ['$T$ = ' num2str(periodVector(5), '%.2f') ' $\mu$s'], ['$T$ = ' num2str(periodVector(6), '%.2f') ' $\mu$s'],...
        ['$T$ = ' num2str(periodVector(7), '%.2f') ' $\mu$s']);
    grid on;
%     xlim([-15 -8]);
    ylim([0 1/sqrt(2)]);
%     formatFig(gcf, [dataPath  'cmin_chirp_dot_' num2str(bandwidthVector(k))], 'en', figProp);
end

save('chirp_dot.mat', 'averageTpr', 'stdTpr', 'averageFpr', 'stdFpr');

rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])
<<<<<<< Updated upstream
rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
    'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
    filesep 'pfa_results' filesep]);
rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
    'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
    filesep]);
=======
if isunix
    rmpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/my_results/pfa_results/');
    rmpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/my_results/');
    
else
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep 'pfa_results' filesep])
    
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep]);
end  

>>>>>>> Stashed changes

%%
%Plot results for chirp for different bandwidths and periods (Pai's
%technique L = 19)

clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '..' filesep '.' filesep 'Misc'])
if isunix
    addpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/pai_results/pfa_results/');
    addpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/pai_results/');
    
else
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep 'pfa_results' filesep])
    
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep]);
end  

load results_det_10.mat;
load resultspfa3.mat;

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
dataPath = ['..' filesep '..' filesep '.' filesep 'figs' filesep '08-12' filesep];

PfaVector = logspace(-12, -2, 41);
bandwidthVector = (2e6:3e6:14e6)/1e6;
periodVector = (8.62e-6:1.48e-6:18.97e-6)*1e6;

JNRVector = -25:0;

for JNRIndex = 1:length(JNRVector)
    for bandwidthIndex = 1:length(bandwidthVector)
        for periodIndex = 1:length(periodVector)
            for thresholdIndex = 1:length(PfaVector)
                x = squeeze(detection_res(bandwidthIndex, periodIndex, JNRIndex, :, :, thresholdIndex));
                tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :) = any(x, 2);
                fn(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :) = ...
                    1 - tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :);
            end
        end
    end
end

tpr = squeeze(tp./(tp+fn));

averageTpr = mean(tpr, 5);
stdTpr = std(tpr, [], 5);

cfun = @(tpr, fpr) sqrt(fpr.^2 + (1-tpr).^2);

for k = 1:length(bandwidthVector)
    for j = 1:length(periodVector)
        for i = 1:length(JNRVector)
            c = cfun(squeeze(averageTpr(k,j,i,:)).', averageFpr(2,:));
            [cMin(k,j,i), idx] = min(c(:));
        end
    end
end

for k = 1:length(bandwidthVector)
    figure;
    
    for j = 1:length(periodVector)
        plot(JNRVector, squeeze(cMin(k,j,:)));
        hold on;
    end
    ylabel('C$_{\mathrm{min}}$');
    xlabel('JNR [dB]');
    legend(['$T$ = ' num2str(periodVector(1), '%.2f') ' $\mu$s'], ['$T$ = ' num2str(periodVector(2), '%.2f') ' $\mu$s'], ['$T$ = ' num2str(periodVector(3), '%.2f') ' $\mu$s'],...
        ['$T$ = ' num2str(periodVector(4), '%.2f') ' $\mu$s'], ['$T$ = ' num2str(periodVector(5), '%.2f') ' $\mu$s'], ['$T$ = ' num2str(periodVector(6), '%.2f') ' $\mu$s'],...
        ['$T$ = ' num2str(periodVector(7), '%.2f') ' $\mu$s']);
    grid on;
%     xlim([-15 -8]);
    ylim([0 1/sqrt(2)]);
%     formatFig(gcf, [dataPath  'cmin_chirp_statistical_19_' num2str(bandwidthVector(k))], 'en', figProp);
end

save('chirp_statistical_19.mat', 'averageTpr', 'stdTpr', 'averageFpr', 'stdFpr');

rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])
if isunix
    rmpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/pai_results/pfa_results/');
    rmpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/pai_results/');
    
else
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep 'pfa_results' filesep])
    
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep]);
end  

%%
%Plot results for chirp for different bandwidths and periods (Pai's
%technique L = 3)

clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '..' filesep '.' filesep 'Misc'])
if isunix
    addpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/pai_results/pfa_results/');
    addpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/pai_results/');
    
else
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep 'pfa_results' filesep])
    
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep]);
end  

load results_det_9.mat;
load resultspfa3.mat;

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
dataPath = ['..' filesep '..' filesep '.' filesep 'figs' filesep '08-12' filesep];

PfaVector = logspace(-12, -2, 41);
bandwidthVector = (2e6:3e6:14e6)/1e6;
periodVector = (8.62e-6:1.48e-6:18.97e-6)*1e6;

JNRVector = -25:0;

for JNRIndex = 1:length(JNRVector)
    for bandwidthIndex = 1:length(bandwidthVector)
        for periodIndex = 1:length(periodVector)
            for thresholdIndex = 1:length(PfaVector)
                x = squeeze(detection_res(bandwidthIndex, periodIndex, JNRIndex, :, thresholdIndex));
                tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :) = squeeze(x);
                fn(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :) = ...
                    1 - tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :);
            end
        end
    end
end

tpr = squeeze(tp./(tp+fn));

averageTpr = mean(tpr, 5);
stdTpr = std(tpr, [], 5);

cfun = @(tpr, fpr) sqrt(fpr.^2 + (1-tpr).^2);

for k = 1:length(bandwidthVector)
    for j = 1:length(periodVector)
        for i = 1:length(JNRVector)
            c = cfun(squeeze(averageTpr(k,j,i,:)).', averageFpr(1,:));
            [cMin(k,j,i), idx] = min(c(:));
        end
    end
end

for k = 1:length(bandwidthVector)
    figure;
    
    for j = 1:length(periodVector)
        plot(JNRVector, squeeze(cMin(k,j,:)));
        hold on;
    end
    ylabel('C$_{\mathrm{min}}$');
    xlabel('JNR [dB]');
    legend(['$T$ = ' num2str(periodVector(1), '%.2f') ' $\mu$s'], ['$T$ = ' num2str(periodVector(2), '%.2f') ' $\mu$s'], ['$T$ = ' num2str(periodVector(3), '%.2f') ' $\mu$s'],...
        ['$T$ = ' num2str(periodVector(4), '%.2f') ' $\mu$s'], ['$T$ = ' num2str(periodVector(5), '%.2f') ' $\mu$s'], ['$T$ = ' num2str(periodVector(6), '%.2f') ' $\mu$s'],...
        ['$T$ = ' num2str(periodVector(7), '%.2f') ' $\mu$s']);
    grid on;
%     xlim([-15 -8]);
    ylim([0 1/sqrt(2)]);
%     formatFig(gcf, [dataPath  'cmin_chirp_statistical_3_' num2str(bandwidthVector(k))], 'en', figProp);
end
save('chirp_statistical_3.mat', 'averageTpr', 'stdTpr', 'averageFpr', 'stdFpr');


rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])
if isunix
    rmpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/pai_results/pfa_results/');
    rmpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/pai_results/');
    
else
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep 'pfa_results' filesep])
    
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep]);
end

%%
%Plot results for chirp for different bandwidths and periods (Pai's
%technique L = 11)

clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '..' filesep '.' filesep 'Misc'])
if isunix
    addpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/pai_results/pfa_results/');
    addpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/pai_results/');
    
else
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep 'pfa_results' filesep])
    
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep]);
end  

load results_det_12.mat;
load resultspfa4.mat;

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
dataPath = ['..' filesep '..' filesep '.' filesep 'figs' filesep '08-12' filesep];

PfaVector = logspace(-12, -2, 41);
bandwidthVector = (2e6:3e6:14e6)/1e6;
periodVector = (8.62e-6:1.48e-6:18.97e-6)*1e6;

JNRVector = -25:0;

for JNRIndex = 1:length(JNRVector)
    for bandwidthIndex = 1:length(bandwidthVector)
        for periodIndex = 1:length(periodVector)
            for thresholdIndex = 1:length(PfaVector)
                x = squeeze(detection_res(bandwidthIndex, periodIndex, JNRIndex, :, thresholdIndex));
                tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :) = squeeze(x);
                fn(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :) = ...
                    1 - tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, :);
            end
        end
    end
end

tpr = squeeze(tp./(tp+fn));

averageTpr = mean(tpr, 5);
stdTpr = std(tpr, [], 5);

cfun = @(tpr, fpr) sqrt(fpr.^2 + (1-tpr).^2);

for k = 1:length(bandwidthVector)
    for j = 1:length(periodVector)
        for i = 1:length(JNRVector)
            c = cfun(squeeze(averageTpr(k,j,i,:)).', averageFpr);
            [cMin(k,j,i), idx] = min(c(:));
        end
    end
end

for k = 1:length(bandwidthVector)
    figure;
    
    for j = 1:length(periodVector)
        plot(JNRVector, squeeze(cMin(k,j,:)));
        hold on;
    end
    ylabel('C$_{\mathrm{min}}$');
    xlabel('JNR [dB]');
    legend(['$T$ = ' num2str(periodVector(1), '%.2f') ' $\mu$s'], ['$T$ = ' num2str(periodVector(2), '%.2f') ' $\mu$s'], ['$T$ = ' num2str(periodVector(3), '%.2f') ' $\mu$s'],...
        ['$T$ = ' num2str(periodVector(4), '%.2f') ' $\mu$s'], ['$T$ = ' num2str(periodVector(5), '%.2f') ' $\mu$s'], ['$T$ = ' num2str(periodVector(6), '%.2f') ' $\mu$s'],...
        ['$T$ = ' num2str(periodVector(7), '%.2f') ' $\mu$s']);
    grid on;
    xlim([-15 -8]);
    ylim([0 1/sqrt(2)]);
    formatFig(gcf, [dataPath  'cmin_chirp_statistical_11_zoom_' num2str(bandwidthVector(k))], 'en', figProp);
end
save('chirp_statistical_11.mat', 'averageTpr', 'stdTpr', 'averageFpr', 'stdFpr');


rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])
if isunix
    rmpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/pai_results/pfa_results/');
    rmpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/pai_results/');
    
else
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep 'pfa_results' filesep])
    
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep]);
end  


%%
%Plot results considering a given period for different bandwidths

clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '..' filesep '.' filesep 'Misc'])


linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
dataPath = ['..' filesep '..' filesep '.' filesep 'figs' filesep '08-12' filesep];

PfaVector = logspace(-12, -2, 41);
bandwidthVector = (2e6:3e6:14e6)/1e6;
periodVector = (8.62e-6:1.48e-6:18.97e-6)*1e6;
JNRVector = -25:0;

cfun = @(tpr, fpr) sqrt(fpr.^2 + (1-tpr).^2);

load chirp_dot.mat;

for k = 1:length(bandwidthVector)
    for j = 1:length(periodVector)
        for i = 1:length(JNRVector)
            c = cfun(squeeze(averageTpr(k,j,i,:)).', averageFpr);
            [cMin(k,j,i), idx] = min(c(:));
        end
    end
end


load chirp_statistical_3.mat;

for k = 1:length(bandwidthVector)
    for j = 1:length(periodVector)
        for i = 1:length(JNRVector)
            c = cfun(squeeze(averageTpr(k,j,i,:)).', averageFpr(1,:));
            [cMin_statistical_3(k,j,i), idx] = min(c(:));
        end
    end
end

load chirp_statistical_19.mat;

for k = 1:length(bandwidthVector)
    for j = 1:length(periodVector)
        for i = 1:length(JNRVector)
            c = cfun(squeeze(averageTpr(k,j,i,:)).', averageFpr(2,:));
            [cMin_statistical_19(k,j,i), idx] = min(c(:));
        end
    end
end

load chirp_statistical_11.mat;

for k = 1:length(bandwidthVector)
    for j = 1:length(periodVector)
        for i = 1:length(JNRVector)
            c = cfun(squeeze(averageTpr(k,j,i,:)).', averageFpr);
            [cMin_statistical_11(k,j,i), idx] = min(c(:));
        end
    end
end


periodIndex = 4;

for k = 1:length(bandwidthVector)
    figure;
    
    plot(JNRVector, squeeze(cMin(k,periodIndex,:)));
    hold on;
    plot(JNRVector, squeeze(cMin_statistical_3(k,periodIndex,:)))
    plot(JNRVector, squeeze(cMin_statistical_11(k,periodIndex,:)))
    plot(JNRVector, squeeze(cMin_statistical_19(k,periodIndex,:)))
    ylabel('C$_{\mathrm{min}}$');
    xlabel('JNR [dB]');
    legend('Dot', 'Statistical, $L = 3$', 'Statistical, $L = 11$', 'Statistical, $L = 19$');
    grid on;
%     xlim([-15 -8]);
    ylim([0 1/sqrt(2)]);
    formatFig(gcf, [dataPath  'cmin_chirp_comparison_' num2str(bandwidthVector(k))], 'en', figProp);
end

rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])