clear;
clc;
close all;

set(groot,'defaultFigureVisible','off')
set(groot,'defaultFigureVisible','on')

addpath(['.' filesep 'data' filesep '04-10']);
addpath(['.' filesep 'Misc'])


linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
expName = 'detect_chirp_sim';
dataPath = ['.' filesep 'figs' filesep '04-10' filesep];

save_fig = false;

threshold = 0.1:0.1:0.9;
load results01.mat;

figure;
for i = 1:3
    
    tp = squeeze(tp_sim_window(:,:,i,1));
    fp = squeeze(fp_sim_window(:,:,i,1));
    fn = squeeze(fn_sim_window(:,:,i,1));
    tn = squeeze(tn_sim_window(:,:,i,1));
    
    fp = fp/3;
    tn = tn/3;
    
    tp = mean(tp, 1);
    fp = mean(fp, 1);
    fn = mean(fn, 1);
    tn = mean(tn, 1);
    
    tpr = tp./(tp+fn);
    
    fpr = fp./(fp+tn);
    
    xAux = linspace(0, 1, numel(fpr));
    yAux = linspace(0, 1, numel(tpr));
    
    plot(fpr, tpr)
    hold on
    
    accuracy1(:,i) = (tp + tn)./(tp + fn + tn + fp);
    [maxAccuracy(i), indexes(i)] = max(accuracy1(:,i), [], 1);

    disp(['TPR: ' num2str(tpr(:, indexes(i))) ' FPR: ' num2str(fpr(:, indexes(i)))])
    disp('------------------------------------')
    disp(['Accuracy: ' num2str(accuracy1(indexes(i),i))])
    disp('------------------------------------')
    disp(['Threshold: ' num2str(threshold(indexes(i)))])
    disp('------------------------------------')
end


plot(xAux, yAux, '--');
box 'off'
grid on
xlabel('Probability of false alarm');
ylabel('Probability of detection');
legend('-10 dB', '-5 dB', '0 dB', 'Random guess', 'Location', 'southeast');
if save_fig
    formatFig(gcf, [dataPath expName '_' 'roc_'  '4'], 'en', figProp); %#ok<*UNRCH>
end

%-------------------------------------------------------

load results02.mat;

figure;
for i = 1:3
    tp = squeeze(tp_sim(:,:,i,1));
    fp = squeeze(fp_sim(:,:,i,1));
    fn = squeeze(fn_sim(:,:,i,1));
    tn = squeeze(tn_sim(:,:,i,1));
    
    fp = fp/2;
    tn = tn/2;
    
    tp = mean(tp, 1);
    fp = mean(fp, 1);
    fn = mean(fn, 1);
    tn = mean(tn, 1);
    
    tpr = tp./(tp+fn);
    
    fpr = fp./(fp+tn);
    
    
    plot(fpr, tpr);
    hold on
   
    accuracy2(:,i) = (tp + tn)./(tp + fn + tn + fp);
    [maxAccuracy(i), indexes(i)] = max(accuracy2(:,i), [], 1);

    disp(['TPR: ' num2str(tpr(:, indexes(i))) ' FPR: ' num2str(fpr(:, indexes(i)))])
    disp('------------------------------------')
    disp(['Accuracy: ' num2str(accuracy2(indexes(i),i))])
    disp('------------------------------------')
    disp(['Threshold: ' num2str(threshold(indexes(i)))])
    disp('------------------------------------')
    
end

xAux = linspace(0, 1, numel(fpr));
yAux = linspace(0, 1, numel(tpr));
plot(xAux, yAux, '--');
box 'off'
grid on
xlabel('Probability of false alarm');
ylabel('Probability of detection');
legend('-10 dB', '-5 dB', '0 dB', 'Random guess', 'Location', 'southeast');

if save_fig
    formatFig(gcf, [dataPath expName '_' 'roc_'  '3'], 'en', figProp); %#ok<*UNRCH>
end

%----------------------------------------------------------


load results03.mat;

figure;
for i = 1:3
    tp = squeeze(tp_window(:,:,i,1));
    fp = squeeze(fp_window(:,:,i,1));
    fn = squeeze(fn_window(:,:,i,1));
    tn = squeeze(tn_window(:,:,i,1));

    fp = fp/3;
    tn = tn/3;

    tp = mean(tp, 1);
    fp = mean(fp, 1);
    fn = mean(fn, 1);
    tn = mean(tn, 1);

    tpr = tp./(tp+fn);

    fpr = fp./(fp+tn);
    
    plot(fpr, tpr);
    hold on
    
    accuracy3(:,i) = (tp + tn)./(tp + fn + tn + fp);
    [maxAccuracy(i), indexes(i)] = max(accuracy3(:,i), [], 1);
    disp(['TPR: ' num2str(tpr(:, indexes(i))) ' FPR: ' num2str(fpr(:, indexes(i)))])
    disp('------------------------------------')
    disp(['Accuracy: ' num2str(accuracy3(indexes(i),i))])
    disp('------------------------------------')
    disp(['Threshold: ' num2str(threshold(indexes(i)))])
    disp('------------------------------------')
    
end

xAux = linspace(0, 1, numel(fpr));
yAux = linspace(0, 1, numel(tpr));
plot(xAux, yAux, '--');
box 'off'
grid on
xlabel('Probability of false alarm');
ylabel('Probability of detection');
legend('-10 dB', '-5 dB', '0 dB', 'Random guess', 'Location', 'southeast');

if save_fig
    formatFig(gcf, [dataPath expName '_' 'roc_'  '2'], 'en', figProp); %#ok<*UNRCH>
end

%------------------------------------------------------


load results04.mat;

figure;
for i = 1:3

    tpAux = squeeze(tp(:,:,i,1));
    fpAux = squeeze(fp(:,:,i,1));
    fnAux = squeeze(fn(:,:,i,1));
    tnAux = squeeze(tn(:,:,i,1));

    fpAux = fpAux/2;
    tnAux = tnAux/2;

    tpAux = mean(tpAux, 1);
    fpAux = mean(fpAux, 1);
    fnAux = mean(fnAux, 1);
    tnAux = mean(tnAux, 1);

    tpr = tpAux./(tpAux+fnAux);

    fpr = fpAux./(fpAux+tnAux);
    
    plot(fpr, tpr)
    hold on

    accuracy4(:,i) = (tpAux + tnAux)./(tpAux + fnAux + tnAux + fpAux);
    [maxAccuracy(i), indexes(i)] = max(accuracy4(:,i), [], 1);

    disp(['TPR: ' num2str(tpr(:, indexes(i))) ' FPR: ' num2str(fpr(:, indexes(i)))])
    disp('------------------------------------')
    disp(['Accuracy: ' num2str(accuracy4(indexes(i),i))])
    disp('------------------------------------')
    disp(['Threshold: ' num2str(threshold(indexes(i)))])
    disp('------------------------------------')

end

xAux = linspace(0, 1, numel(fpr));
yAux = linspace(0, 1, numel(tpr));
plot(xAux, yAux, '--');
box 'off'
grid on
xlabel('Probability of false alarm');
ylabel('Probability of detection');
legend('-10 dB', '-5 dB', '0 dB', 'Random guess', 'Location', 'southeast');

if save_fig
    formatFig(gcf, [dataPath expName '_' 'roc_'  '1'], 'en', figProp); %#ok<*UNRCH>
end
    
    
rmpath(['.' filesep 'Misc'])
rmpath(['.' filesep 'data' filesep '04-10']);