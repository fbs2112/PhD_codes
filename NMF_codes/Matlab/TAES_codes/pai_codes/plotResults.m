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
    fp = zeros(length(PfaVector), monteCarloLoops);
    tn = zeros(length(PfaVector), monteCarloLoops);
    for j = 1:length(PfaVector)
        for i = 1:monteCarloLoops
       
        fp(j,i) = sum(x(i,:,j));
        tn(j,i) = size(x, 2) - fp(j,i);
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
%%
%ROC curve under -17 dB JNR for cw interference
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

load results05.mat;
load pfa_data.mat;

PfaVector = logspace(-5, 0, 17);
monteCarloLoops = 100;
pdIndex = 26;

JNRIdx = 7;

for j = 1:length(PfaVector)
    for i = 1:monteCarloLoops
        
        tp(j,i) = detection_res(JNRIdx,i,pdIndex,j);
        fn(j,i) = 1 - tp(j,i);
    end
end
tpr = tp./(tp+fn);

averageTprPai = squeeze(mean(tpr, 2));
stdTpr = squeeze(std(tpr, [], 2));

averageFprPai = averageFpr;

figure;
plot(averageFprPai(2,:) , averageTprPai.')

ylabel('Probability of false alarm');
xlabel('$\bar{\gamma}$');
xlim([min(PfaVector) max(PfaVector)]);
grid on;


figure;
loglog(PfaVector , averageTprPai.')

ylabel('Probability of false alarm');
xlabel('$\bar{\gamma}$');
xlim([min(PfaVector) max(PfaVector)]);
ylim([1e-5 1e-0]);
grid on;

figure;
semilogx(averageFprPai(2,:) , averageTprPai.')
grid on;
hold on;

rmpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'pai_results']);
%------------------------------------------------------------------------------------------------

addpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results']);  

load pfa_data;

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);

monteCarloLoops = 100;

thresholdVector = 0.1:0.05:0.9;
window_median_length_vector = 51:50:401;
periodVector = 0;
bandwidthVector = 0;
JNRVector = 1;

tp = zeros(length(bandwidthVector), length(periodVector), length(JNRVector), length(thresholdVector), length(window_median_length_vector), monteCarloLoops);
fn = zeros(length(bandwidthVector), length(periodVector), length(JNRVector), length(thresholdVector), length(window_median_length_vector), monteCarloLoops);

for JNRIndex = JNRIdx:JNRIdx%length(JNRVector)
    load(['results05_' num2str(JNRIndex)]);
    for bandwidthIndex = 1:length(bandwidthVector)
        for periodIndex = 1:length(periodVector)
            for thresholdIndex = 1:length(thresholdVector)
                for window_median_length_index = 1:length(window_median_length_vector)
                    x = squeeze(detection_res(:, bandwidthIndex, periodIndex, 1, thresholdIndex, window_median_length_index,:));
                    tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, window_median_length_index, :) = sum(x, 2);
                    fn(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, window_median_length_index, :) = ...
                        size(x, 2) - tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, window_median_length_index, :);
                end
            end
        end
    end
end

tpr = squeeze(tp./(tp+fn));

averageTpr = mean(tpr, 4);
stdTpr = std(tpr, [], 4);

for i = 1:size(averageTpr, 3)
    semilogx(averageFpr(i,:), squeeze(averageTpr(JNRIdx,:,i)));
end

xlabel('Probability of false alarm');
ylabel('Probability of detection');
xlim([min(PfaVector) max(PfaVector)]);
ylim([1e-2 1e-0]);

rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])
rmpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results']);  

%%
%ROC curve under -17 dB JNR for cw interference with silence periods
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

load results06.mat;
load pfa_data.mat;

PfaVector = logspace(-5, 0, 17);
monteCarloLoops = 100;
pdIndex = 36;

JNRIdx = 6;

for j = 1:length(PfaVector)
    for i = 1:monteCarloLoops
        
        tp(j,i) = detection_res(JNRIdx,i,pdIndex,j);
        fn(j,i) = 1 - tp(j,i);
    end
end
tpr = tp./(tp+fn);

averageTprPai = squeeze(mean(tpr, 2));
stdTpr = squeeze(std(tpr, [], 2));

averageFprPai = averageFpr;

figure;
plot(averageFprPai(2,:) , averageTprPai.')

ylabel('Probability of false alarm');
xlabel('$\bar{\gamma}$');
xlim([min(PfaVector) max(PfaVector)]);
grid on;


figure;
loglog(PfaVector , averageTprPai.')

ylabel('Probability of false alarm');
xlabel('$\bar{\gamma}$');
xlim([min(PfaVector) max(PfaVector)]);
ylim([1e-5 1e-0]);
grid on;

figure;
semilogx(averageFprPai(2,:) , averageTprPai.')
grid on;
hold on;

rmpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'pai_results']);
%------------------------------------------------------------------------------------------------

addpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results']);  

load pfa_data;

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);

monteCarloLoops = 100;

thresholdVector = 0.1:0.05:0.9;
window_median_length_vector = 51:50:401;
periodVector = 0;
bandwidthVector = 0;
JNRVector = 1;

window_length = 98;
onset = 528;
offset = 4623;

tp = zeros(length(bandwidthVector), length(periodVector), length(JNRVector), length(thresholdVector), length(window_median_length_vector), monteCarloLoops);
fn = zeros(length(bandwidthVector), length(periodVector), length(JNRVector), length(thresholdVector), length(window_median_length_vector), monteCarloLoops);

for JNRIndex = JNRIdx:JNRIdx%length(JNRVector)
    load(['results06_' num2str(JNRIndex)]);
    for bandwidthIndex = 1:length(bandwidthVector)
        for periodIndex = 1:length(periodVector)
            for thresholdIndex = 1:length(thresholdVector)
                for window_median_length_index = 1:length(window_median_length_vector)
                    x = squeeze(detection_res(:, bandwidthIndex, periodIndex, 1, thresholdIndex, window_median_length_index,onset  + round(window_length*2):offset - round(window_length*2)));
                    tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, window_median_length_index, :) = sum(x, 2);
                    fn(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, window_median_length_index, :) = ...
                        size(x, 2) - tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, window_median_length_index, :);
                end
            end
        end
    end
end

tpr = squeeze(tp./(tp+fn));

averageTpr = mean(tpr, 4);
stdTpr = std(tpr, [], 4);

for i = 1:size(averageTpr, 3)
    semilogx(averageFpr(i,:), squeeze(averageTpr(JNRIdx,:,i)));
end

xlabel('Probability of false alarm');
ylabel('Probability of detection');
xlim([min(PfaVector) max(PfaVector)]);
ylim([1e-2 1e-0]);

rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])
rmpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results']);  
%%
%ROC curve under -17 dB JNR for chirp interference
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

load results07.mat;
load pfa_data.mat;

PfaVector = logspace(-5, 0, 17);
monteCarloLoops = 100;
pdIndex = 300;
JNRIdx = 10;

for j = 1:length(PfaVector)
    for i = 1:monteCarloLoops
        
        tp(j,i) = detection_res(JNRIdx,i,pdIndex,j);
        fn(j,i) = 1 - tp(j,i);
    end
end
tpr = tp./(tp+fn);

averageTprPai = squeeze(mean(tpr, 2));
stdTpr = squeeze(std(tpr, [], 2));

averageFprPai = averageFpr;

figure;
plot(averageFprPai(2,:) , averageTprPai.')

ylabel('Probability of false alarm');
xlabel('$\bar{\gamma}$');
xlim([min(PfaVector) max(PfaVector)]);
grid on;

figure;
loglog(PfaVector , averageTprPai.')

ylabel('Probability of false alarm');
xlabel('$\bar{\gamma}$');
xlim([min(PfaVector) max(PfaVector)]);
ylim([1e-5 1e-0]);
grid on;

figure;
semilogx(averageFprPai(1,:) , averageTprPai.')
grid on;
hold on;

rmpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'pai_results']);
%------------------------------------------------------------------------------------------------

addpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results']);  

load pfa_data;

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);

monteCarloLoops = 100;

thresholdVector = 0.1:0.05:0.9;
window_median_length_vector = 51:50:401;
periodVector = 0;
bandwidthVector = 0;
JNRVector = 1;

tp = zeros(length(bandwidthVector), length(periodVector), length(JNRVector), length(thresholdVector), length(window_median_length_vector), monteCarloLoops);
fn = zeros(length(bandwidthVector), length(periodVector), length(JNRVector), length(thresholdVector), length(window_median_length_vector), monteCarloLoops);

for JNRIndex = JNRIdx:JNRIdx%length(JNRVector)
    load(['results07_' num2str(JNRIndex)]);
    for bandwidthIndex = 1:length(bandwidthVector)
        for periodIndex = 1:length(periodVector)
            for thresholdIndex = 1:length(thresholdVector)
                for window_median_length_index = 1:length(window_median_length_vector)
                    x = squeeze(detection_res(:, bandwidthIndex, periodIndex, 1, thresholdIndex, window_median_length_index,:));
                    tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, window_median_length_index, :) = sum(x, 2);
                    fn(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, window_median_length_index, :) = ...
                        size(x, 2) - tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, window_median_length_index, :);
                end
            end
        end
    end
end

tpr = squeeze(tp./(tp+fn));

averageTpr = mean(tpr, 4);
stdTpr = std(tpr, [], 4);

for i = 1:size(averageTpr, 3)
    semilogx(averageFpr(i,:), squeeze(averageTpr(JNRIdx,:,i)));
end

xlabel('Probability of false alarm');
ylabel('Probability of detection');
xlim([min(PfaVector) max(PfaVector)]);
ylim([1e-2 1e-0]);

rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])
rmpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results']);  

%%
%ROC curve under -17 dB JNR for chirp interference with silence
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

load results08.mat;
load pfa_data.mat;

PfaVector = logspace(-5, 0, 17);
monteCarloLoops = 100;
pdIndex = 400;
JNRIdx = 10;

for j = 1:length(PfaVector)
    for i = 1:monteCarloLoops
        
        tp(j,i) = detection_res(JNRIdx,i,pdIndex,j);
        fn(j,i) = 1 - tp(j,i);
    end
end
tpr = tp./(tp+fn);

averageTprPai = squeeze(mean(tpr, 2));
stdTpr = squeeze(std(tpr, [], 2));

averageFprPai = averageFpr;

figure;
plot(averageFprPai(2,:) , averageTprPai.')

ylabel('Probability of false alarm');
xlabel('$\bar{\gamma}$');
xlim([min(PfaVector) max(PfaVector)]);
grid on;

figure;
loglog(PfaVector , averageTprPai.')

ylabel('Probability of false alarm');
xlabel('$\bar{\gamma}$');
xlim([min(PfaVector) max(PfaVector)]);
ylim([1e-5 1e-0]);
grid on;

figure;
semilogx(averageFprPai(1,:) , averageTprPai.')
grid on;
hold on;

rmpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'pai_results']);
%------------------------------------------------------------------------------------------------

addpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results']);  

load pfa_data;

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);

window_length = 98;
onset = 528;
offset = 4623;

monteCarloLoops = 100;

thresholdVector = 0.1:0.05:0.9;
window_median_length_vector = 51:50:401;
periodVector = 0;
bandwidthVector = 0;
JNRVector = 1;

tp = zeros(length(bandwidthVector), length(periodVector), length(JNRVector), length(thresholdVector), length(window_median_length_vector), monteCarloLoops);
fn = zeros(length(bandwidthVector), length(periodVector), length(JNRVector), length(thresholdVector), length(window_median_length_vector), monteCarloLoops);

for JNRIndex = JNRIdx:JNRIdx%length(JNRVector)
    load(['results08_' num2str(JNRIndex)]);
    for bandwidthIndex = 1:length(bandwidthVector)
        for periodIndex = 1:length(periodVector)
            for thresholdIndex = 1:length(thresholdVector)
                for window_median_length_index = 1:length(window_median_length_vector)
                    x = squeeze(detection_res(:, bandwidthIndex, periodIndex, 1, thresholdIndex, window_median_length_index,onset  + round(window_length*2):offset - round(window_length*2)));
                    tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, window_median_length_index, :) = sum(x, 2);
                    fn(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, window_median_length_index, :) = ...
                        size(x, 2) - tp(bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, window_median_length_index, :);
                end
            end
        end
    end
end

tpr = squeeze(tp./(tp+fn));

averageTpr = mean(tpr, 4);
stdTpr = std(tpr, [], 4);

for i = 1:size(averageTpr, 3)
    semilogx(averageFpr(i,:), squeeze(averageTpr(JNRIdx,:,i)));
end

xlabel('Probability of false alarm');
ylabel('Probability of detection');
xlim([min(PfaVector) max(PfaVector)]);
ylim([1e-2 1e-0]);

rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])
rmpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results']);  
