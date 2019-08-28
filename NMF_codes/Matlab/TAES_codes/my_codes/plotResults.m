clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '..' filesep '.' filesep 'Misc'])
addpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results']);  

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);

load results17.mat;

thresholdVector = 0.1:0.05:0.9;
window_median_length_vector = 51:50:401;
window_median_length_vector = 0;


for k = 1:length(thresholdVector)
    for j = 1:length(window_median_length_vector)
       
        x = squeeze(detection_res(:,k));
        fp(:,k) = x;
        tn(:,k) = 1 - fp(:,k) ;
    end
end

fpr = fp./(fp+tn);

averageFpr = squeeze(mean(fpr, 1));
stdFpr = squeeze(std(fpr, [], 1));

figure;
% for i = 1:size(averageFpr, 1)
plot(thresholdVector, averageFpr)
%     hold on;
% end

ylabel('Probability of false alarm');
xlabel('$\bar{\gamma}$');
xlim([min(thresholdVector) max(thresholdVector)]);

% legend('$L_{\mathrm{med}} = 51$', '$L_{\mathrm{med}} = 101$', '$L_{\mathrm{med}} = 151$','$L_{\mathrm{med}} = 201$',...
%     '$L_{\mathrm{med}} = 251$', '$L_{\mathrm{med}} = 301$', '$L_{\mathrm{med}} = 351$', '$L_{\mathrm{med}} = 401$');
grid on;

figure;
% for i = 1:size(averageFpr, 1)
    loglog(thresholdVector, averageFpr)
%     hold on;
% end

ylabel('Probability of false alarm');
xlabel('$\bar{\gamma}$');
xlim([min(thresholdVector) max(thresholdVector)]);

% legend('$L_{\mathrm{med}} = 51$', '$L_{\mathrm{med}} = 101$', '$L_{\mathrm{med}} = 151$','$L_{\mathrm{med}} = 201$',...
%     '$L_{\mathrm{med}} = 251$', '$L_{\mathrm{med}} = 301$', '$L_{\mathrm{med}} = 351$', '$L_{\mathrm{med}} = 401$');
grid on;

save(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results' filesep 'pfa_data_median_full_128.mat'], 'averageFpr', 'stdFpr')

rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])
rmpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results']);  

%%
clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '..' filesep '.' filesep 'Misc'])
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
JNRVector = -20:0;

for JNRIndex = 1:length(JNRVector)
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

figure;
for i = 1:size(averageTpr, 3)
    plot(thresholdVector, averageTpr(6,:,i))
    hold on
end

ylabel('Probability of detection');
xlabel('$\bar{\gamma}$');
xlim([min(thresholdVector) max(thresholdVector)]);
legend('$L_{\mathrm{med}} = 51$', '$L_{\mathrm{med}} = 101$', '$L_{\mathrm{med}} = 151$','$L_{\mathrm{med}} = 201$',...
    '$L_{\mathrm{med}} = 251$', '$L_{\mathrm{med}} = 301$', '$L_{\mathrm{med}} = 351$', '$L_{\mathrm{med}} = 401$');
grid on;

figure;
for i = 1:size(averageTpr, 3)
    loglog(thresholdVector, averageTpr(6,:,i))
    hold on
end

ylabel('Probability of detection');
xlabel('$\bar{\gamma}$');
xlim([min(thresholdVector) max(thresholdVector)]);
legend('$L_{\mathrm{med}} = 51$', '$L_{\mathrm{med}} = 101$', '$L_{\mathrm{med}} = 151$','$L_{\mathrm{med}} = 201$',...
    '$L_{\mathrm{med}} = 251$', '$L_{\mathrm{med}} = 301$', '$L_{\mathrm{med}} = 351$', '$L_{\mathrm{med}} = 401$');
grid on;
    
cfun = @(tpr, fpr) sqrt(fpr.^2 + (1-tpr).^2);

for i = 1:length(JNRVector)
    c = cfun(squeeze(averageTpr(i,:,:)).', averageFpr);
    cMin(i) = min(c(:));
end

figure;
plot(JNRVector, cMin)

%Plot ROC
figure;
for i = 1:size(averageTpr, 3)
    loglog(averageFpr(i,:), averageTpr(4, :, i));
    hold on;
end
legend('$L_{\mathrm{med}} = 51$', '$L_{\mathrm{med}} = 101$', '$L_{\mathrm{med}} = 151$','$L_{\mathrm{med}} = 201$',...
    '$L_{\mathrm{med}} = 251$', '$L_{\mathrm{med}} = 301$', '$L_{\mathrm{med}} = 351$', '$L_{\mathrm{med}} = 401$');
xlim([1e-5 1e0]);
ylim([1e-5 1e0]);
grid on;

rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])
rmpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results']); 

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
addpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results']);  

load pfa_data_median_full_256;

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);

load results20.mat;

monteCarloLoops = 100;

thresholdVector = 0.1:0.05:0.9;
window_median_length_vector = 0;
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

JNRIdx = 1;
figure;
plot(thresholdVector, averageTpr(JNRIdx,:))

ylabel('Probability of detection');
xlabel('$\bar{\gamma}$');
xlim([min(thresholdVector) max(thresholdVector)]);
grid on;

figure;
loglog(thresholdVector, averageTpr(JNRIdx,:))

ylabel('Probability of detection');
xlabel('$\bar{\gamma}$');
xlim([min(thresholdVector) max(thresholdVector)]);
grid on;
    
% for i = 1:size(averageTpr, 1)
%     figure;
%     plot(averageFpr, averageTpr(i,:));
%     hold on 
%     plot(linspace(0, 1, numel(averageFpr)), linspace(0, 1, numel(averageFpr)));
%     title(['JNR: ' num2str(JNRVector(i))]);
%     ylabel('Probability of detection');
%     xlabel('Probability of false alarm');
%     grid on;
% end
cfun = @(tpr, fpr) sqrt(fpr.^2 + (1-tpr).^2);

for i = 1:length(JNRVector)
    c = cfun(averageTpr(i,:), averageFpr);
    cMin(i) = min(c(:));
end

figure;
plot(JNRVector, cMin)
grid on;
ylabel('C$_{\mathrm{min}}$');
xlabel('JNR [dB]');
rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])
rmpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results']); 

%%
%Plot results for chirp considering the new approach of using a large median
%window

clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '..' filesep '.' filesep 'Misc'])
addpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results']);  

load pfa_data_median_full_64_lowT.mat;

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);

load results48.mat;

monteCarloLoops = 100;

thresholdVector = -0.1:0.01:0.2;
window_median_length_vector = 0;
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

JNRIdx = 1;
% figure;
% plot(thresholdVector, averageTpr(JNRIdx,:))
% 
% ylabel('Probability of detection');
% xlabel('$\bar{\gamma}$');
% xlim([min(thresholdVector) max(thresholdVector)]);
% grid on;
% 
% figure;
% loglog(thresholdVector, averageTpr(JNRIdx,:))
% 
% ylabel('Probability of detection');
% xlabel('$\bar{\gamma}$');
% xlim([min(thresholdVector) max(thresholdVector)]);
% grid on;
%     
% for i = 1:size(averageTpr, 1)
%     figure;
%     plot(averageFpr, averageTpr(i,:));
%     hold on 
%     plot(linspace(0, 1, numel(averageFpr)), linspace(0, 1, numel(averageFpr)));
%     title(['JNR: ' num2str(JNRVector(i))]);
%     ylabel('Probability of detection');
%     xlabel('Probability of false alarm');
%     grid on;
% end

cfun = @(tpr, fpr) sqrt(fpr.^2 + (1-tpr).^2);

for i = 1:length(JNRVector)
    c = cfun(averageTpr(i,:), averageFpr);
    cMin(i) = min(c(:));
end

figure;
plot(JNRVector, cMin)
grid on;
ylabel('C$_{\mathrm{min}}$');
xlabel('JNR [dB]');
ylim([0 1/sqrt(2)])
xlim([min(JNRVector) max(JNRVector)]);
rmpath(['..' filesep '..' filesep '.' filesep 'Misc'])
rmpath(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results']); 