clear;
clc;
close all;

addpath(['..' filesep '.' filesep 'Misc'])


linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
expName = 'detect_chirp_sim_gaussian';
dataPath = ['..' filesep '.' filesep 'figs' filesep '05-02' filesep];

save_fig = false;

threshold = 0.1:0.1:0.9;
stdVector = 7:2:13;

addpath(['..' filesep '.' filesep 'data' filesep '05-02']);    

load results01.mat;

for j = 1:4
    disp(['sigma: ' num2str(stdVector(j))])
    figure;
    for i = 1:5
        tp = squeeze(tp_sim(:,:,i,j));
        fp = squeeze(fp_sim(:,:,i,j));
        fn = squeeze(fn_sim(:,:,i,j));
        tn = squeeze(tn_sim(:,:,i,j));

        fp = fp/2;
        tn = tn/2;

        
        accuracyStd = std((tp + tn)./(tp + fn + tn + fp), 1);
        
        tpr = tp./(tp+fn);
        fpr = fp./(fp+tn);
        
        tprStd = std(tpr, 1);
        fprStd = std(fpr, 1);
        
        tp = mean(tp, 1);
        fp = mean(fp, 1);
        fn = mean(fn, 1);
        tn = mean(tn, 1);

        tpr = mean(tpr, 1);
        fpr = mean(fpr, 1);

        plot(fpr, tpr);
        hold on

        accuracy2(:,i) = (tp + tn)./(tp + fn + tn + fp);
        [maxAccuracy(i), indexes(i)] = max(accuracy2(:,i), [], 1);
        
        disp('------------------------------------')
        disp(['TPR: ' num2str(tpr(:, indexes(i))) '+' num2str(tprStd(:, indexes(i))) ' FPR: ' num2str(fpr(:, indexes(i))) '+' num2str(fprStd(:, indexes(i)))])
        disp('------------------------------------')
        disp(['Accuracy: ' num2str(accuracy2(indexes(i),i)) '+' num2str(accuracyStd(indexes(i)))])
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
    legend('-20 dB', '-15 dB', '-10 dB', '-5 dB', '0 dB', 'Random guess', 'Location', 'southeast');

    if save_fig
        formatFig(gcf, [dataPath expName '_' 'roc_'  '1_' num2str(j)], 'en', figProp); %#ok<*UNRCH>
    end
%-------------------------------------------------------

end

clc;
load results02.mat;

for j = 1:4
    
    
    figure;
    for i = 1:5
        tp = squeeze(tp_sim_window(:,:,i,j));
        fp = squeeze(fp_sim_window(:,:,i,j));
        fn = squeeze(fn_sim_window(:,:,i,j));
        tn = squeeze(tn_sim_window(:,:,i,j));
        
        fp = fp/3;
        tn = tn/3;
        
        accuracyStd = std((tp + tn)./(tp + fn + tn + fp), 1);
        
        tpr = tp./(tp+fn);
        fpr = fp./(fp+tn);
        
        tprStd = std(tpr, 1);
        fprStd = std(fpr, 1);
        
        tp = mean(tp, 1);
        fp = mean(fp, 1);
        fn = mean(fn, 1);
        tn = mean(tn, 1);

        tpr = mean(tpr, 1);
        fpr = mean(fpr, 1);

        plot(fpr, tpr);
        hold on
        
        accuracy2(:,i) = (tp + tn)./(tp + fn + tn + fp);
        [maxAccuracy(i), indexes(i)] = max(accuracy2(:,i), [], 1);
        
        disp('------------------------------------')
        disp(['TPR: ' num2str(tpr(:, indexes(i))) '+' num2str(tprStd(:, indexes(i))) ' FPR: ' num2str(fpr(:, indexes(i))) '+' num2str(fprStd(:, indexes(i)))])
        disp('------------------------------------')
        disp(['Accuracy: ' num2str(accuracy2(indexes(i),i)) '+' num2str(accuracyStd(indexes(i)))])
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
    legend('-20 dB', '-15 dB', '-10 dB', '-5 dB', '0 dB', 'Random guess', 'Location', 'southeast');
    
    if save_fig
        formatFig(gcf, [dataPath expName '_' 'roc_'  '2_' num2str(j)], 'en', figProp); %#ok<*UNRCH>
    end
    
end
%-------------------------------------------------------
clc;

load results03.mat;

expName = 'detect_chirp_sim_inner_median';

figure;
for i = 1:5
    tp = squeeze(tp_sim(:,:,i,1));
    fp = squeeze(fp_sim(:,:,i,1));
    fn = squeeze(fn_sim(:,:,i,1));
    tn = squeeze(tn_sim(:,:,i,1));
    
    fp = fp/2;
    tn = tn/2;
   
    accuracyStd = std((tp + tn)./(tp + fn + tn + fp), 1);
    
    tpr = tp./(tp+fn);
    fpr = fp./(fp+tn);
    
    tprStd = std(tpr, 1);
    fprStd = std(fpr, 1);
    
    tp = mean(tp, 1);
    fp = mean(fp, 1);
    fn = mean(fn, 1);
    tn = mean(tn, 1);
    
    tpr = mean(tpr, 1);
    fpr = mean(fpr, 1);
    
    plot(fpr, tpr);
    hold on
    
    accuracy2(:,i) = (tp + tn)./(tp + fn + tn + fp);
    [maxAccuracy(i), indexes(i)] = max(accuracy2(:,i), [], 1);
    
    disp('------------------------------------')
    disp(['TPR: ' num2str(tpr(:, indexes(i))) '+' num2str(tprStd(:, indexes(i))) ' FPR: ' num2str(fpr(:, indexes(i))) '+' num2str(fprStd(:, indexes(i)))])
    disp('------------------------------------')
    disp(['Accuracy: ' num2str(accuracy2(indexes(i),i)) '+' num2str(accuracyStd(indexes(i)))])
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
legend('-20 dB', '-15 dB', '-10 dB', '-5 dB', '0 dB', 'Random guess', 'Location', 'southeast');

if save_fig
    formatFig(gcf, [dataPath expName '_' 'roc_'  '1_'], 'en', figProp); %#ok<*UNRCH>
end
%-------------------------------------------------------

clc;

load results04.mat;

expName = 'detect_chirp_sim_inner_median_window';

figure;
for i = 1:5
    tp = squeeze(tp_sim_window(:,:,i,1));
    fp = squeeze(fp_sim_window(:,:,i,1));
    fn = squeeze(fn_sim_window(:,:,i,1));
    tn = squeeze(tn_sim_window(:,:,i,1));
    
    fp = fp/3;
    tn = tn/3;
    
    
    accuracyStd = std((tp + tn)./(tp + fn + tn + fp), 1);
    
    tpr = tp./(tp+fn);
    fpr = fp./(fp+tn);
    
    tprStd = std(tpr, 1);
    fprStd = std(fpr, 1);
    
    tp = mean(tp, 1);
    fp = mean(fp, 1);
    fn = mean(fn, 1);
    tn = mean(tn, 1);
    
    tpr = mean(tpr, 1);
    fpr = mean(fpr, 1);
    
    plot(fpr, tpr);
    hold on
    
    accuracy2(:,i) = (tp + tn)./(tp + fn + tn + fp);
    [maxAccuracy(i), indexes(i)] = max(accuracy2(:,i), [], 1);
    
    disp('------------------------------------')
    disp(['TPR: ' num2str(tpr(:, indexes(i))) '+' num2str(tprStd(:, indexes(i))) ' FPR: ' num2str(fpr(:, indexes(i))) '+' num2str(fprStd(:, indexes(i)))])
    disp('------------------------------------')
    disp(['Accuracy: ' num2str(accuracy2(indexes(i),i)) '+' num2str(accuracyStd(indexes(i)))])
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
legend('-20 dB', '-15 dB', '-10 dB', '-5 dB', '0 dB', 'Random guess', 'Location', 'southeast');

if save_fig
    formatFig(gcf, [dataPath expName '_' 'roc_'  '1_'], 'en', figProp); %#ok<*UNRCH>
end
%-------------------------------------------------------
clc;

load results07.mat;
expName = 'detect_chirp_sim_gaussian_median';


for j = 1:4
    disp(['sigma: ' num2str(stdVector(j))])
    figure;
    for i = 1:5
        tp = squeeze(tp_sim(:,:,i,j));
        fp = squeeze(fp_sim(:,:,i,j));
        fn = squeeze(fn_sim(:,:,i,j));
        tn = squeeze(tn_sim(:,:,i,j));

        fp = fp/2;
        tn = tn/2;
        
        accuracyStd = std((tp + tn)./(tp + fn + tn + fp), 1);
        
        tpr = tp./(tp+fn);
        fpr = fp./(fp+tn);
        
        tprStd = std(tpr, 1);
        fprStd = std(fpr, 1);
        
        tp = mean(tp, 1);
        fp = mean(fp, 1);
        fn = mean(fn, 1);
        tn = mean(tn, 1);

        tpr = mean(tpr, 1);
        fpr = mean(fpr, 1);

        plot(fpr, tpr);
        hold on

        accuracy2(:,i) = (tp + tn)./(tp + fn + tn + fp);
        [maxAccuracy(i), indexes(i)] = max(accuracy2(:,i), [], 1);
        
        disp('------------------------------------')
        disp(['TPR: ' num2str(tpr(:, indexes(i))) '+' num2str(tprStd(:, indexes(i))) ' FPR: ' num2str(fpr(:, indexes(i))) '+' num2str(fprStd(:, indexes(i)))])
        disp('------------------------------------')
        disp(['Accuracy: ' num2str(accuracy2(indexes(i),i)) '+' num2str(accuracyStd(indexes(i)))])
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
    legend('-20 dB', '-15 dB', '-10 dB', '-5 dB', '0 dB', 'Random guess', 'Location', 'southeast');

    if save_fig
        formatFig(gcf, [dataPath expName '_' 'roc_'  '1_' num2str(j)], 'en', figProp); %#ok<*UNRCH>
    end
%-------------------------------------------------------

end

clc;

load results06.mat;
expName = 'detect_chirp_sim_gaussian_window_median';
save_fig = true;
for j = 1:4
    
    
    figure;
    for i = 1:5
        tp = squeeze(tp_sim_window(:,:,i,j));
        fp = squeeze(fp_sim_window(:,:,i,j));
        fn = squeeze(fn_sim_window(:,:,i,j));
        tn = squeeze(tn_sim_window(:,:,i,j));
        
        fp = fp/3;
        tn = tn/3;
        
        accuracyStd = std((tp + tn)./(tp + fn + tn + fp), 1);
        
        tpr = tp./(tp+fn);
        fpr = fp./(fp+tn);
        
        tprStd = std(tpr, 1);
        fprStd = std(fpr, 1);
        
        tp = mean(tp, 1);
        fp = mean(fp, 1);
        fn = mean(fn, 1);
        tn = mean(tn, 1);

        tpr = mean(tpr, 1);
        fpr = mean(fpr, 1);

        plot(fpr, tpr);
        hold on
        
        accuracy2(:,i) = (tp + tn)./(tp + fn + tn + fp);
        [maxAccuracy(i), indexes(i)] = max(accuracy2(:,i), [], 1);
        
        disp('------------------------------------')
        disp(['TPR: ' num2str(tpr(:, indexes(i))) '+' num2str(tprStd(:, indexes(i))) ' FPR: ' num2str(fpr(:, indexes(i))) '+' num2str(fprStd(:, indexes(i)))])
        disp('------------------------------------')
        disp(['Accuracy: ' num2str(accuracy2(indexes(i),i)) '+' num2str(accuracyStd(indexes(i)))])
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
    legend('-20 dB', '-15 dB', '-10 dB', '-5 dB', '0 dB', 'Random guess', 'Location', 'southeast');
    
    if save_fig
        formatFig(gcf, [dataPath expName '_' 'roc_'  '2_' num2str(j)], 'en', figProp); %#ok<*UNRCH>
    end
    
end
%-------------------------------------------------------


rmpath(['..' filesep '.' filesep 'data' filesep '05-02']);    
rmpath(['..' filesep '.' filesep 'Misc'])
