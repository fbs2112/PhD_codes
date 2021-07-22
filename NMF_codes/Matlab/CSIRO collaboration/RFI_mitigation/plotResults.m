clear;
clc;
close all;
addpath(['..' filesep '..' filesep 'Misc'])

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

linewidth = 1.5;
fontname = 'Times';
fontsize = 26;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
dataPath = ['..' filesep '..' filesep 'figs' filesep '2021-05-20' filesep];

frequency = '1152';
fileName = '2018-09-01-09_52_55_0000000000000000.000000';
resultsFileNumber = 'results3';
year = '2018';
filePath2 = ['E:' filesep 'Parkes data' filesep year filesep frequency filesep 'signal frames' filesep];

matObj = matfile(['..' filesep 'Labels' filesep year filesep frequency filesep fileName '_labels.mat']);
trueLabels = find(matObj.interferenceDetFlag == 1 | matObj.interferenceDetFlag == 3 | matObj.interferenceDetFlag == 5 | matObj.interferenceDetFlag == 6);

% trueLabels = find(interferenceDetFlag);
% trueLabels = sort(trueLabels);
% trueLabels = trueLabels(trueLabels < 100);
numberOfSignalFrames = 76;

Fs = 128e6;
nfft = 256;
nperseg = 256;
overlap = nperseg/2;
window = ones(nperseg, 1);
mindB = 0;
maxdB = 70;

for loopIndex = 21:numberOfSignalFrames
    disp(['Loop index: ' num2str(loopIndex)]);
    load([filePath2 fileName '_' num2str(trueLabels(loopIndex)) '.mat']);
    matObjRes = matfile(['.' filesep 'data' filesep year filesep frequency filesep resultsFileNumber '.mat']);
    
    [PxxAOr, ~, ~] = spectrogram(parkesSignalA, window, overlap, nfft, Fs, 'centered');
    [PxxBOr, ~, ~] = spectrogram(parkesSignalB, window, overlap, nfft, Fs, 'centered');
    
    [PxxA, ~, ~] = spectrogram(squeeze(matObjRes.xHatA(:,2,1,trueLabels(loopIndex))), window, overlap, nfft, Fs, 'centered');
    [PxxB, f, t] = spectrogram(squeeze(matObjRes.xHatB(:,2,1,trueLabels(loopIndex))), window, overlap, nfft, Fs, 'centered');
    
%     [PxxARFI, ~, ~] = spectrogram(squeeze(matObjRes.xHatA(:,1,1,trueLabels(loopIndex))), window, overlap, nfft, Fs, 'centered');
%     [PxxBRFI, f, t] = spectrogram(squeeze(matObjRes.xHatB(:,1,1,trueLabels(loopIndex))), window, overlap, nfft, Fs, 'centered');
    
    figure;
    surf(t*1e3, f/1e6 + str2double(frequency), 20*log10(abs(PxxAOr)), 'EdgeColor', 'none');
    axis xy;
    view(0, 90);
    xlim([min(t) max(t)]*1e3);
    ylim([min(f) max(f)]/1e6 + str2double(frequency));
    xlabel('Time [ms]');
    ylabel('Frequency [MHz]');
    ax = gca;
    set(ax, 'clim', [mindB maxdB], 'colormap', turbo);
    c = colorbar;
    c.Label.String = '[dB]';
    c.TickLabelInterpreter = 'latex';
    
    formatFig(gcf, [dataPath 'original_' num2str(loopIndex)], 'en', figProp);
    formatFig(gcf, [dataPath 'original_' num2str(loopIndex)], 'en', figProp);
    
    figure;
    surf(t*1e3, f/1e6 + str2double(frequency), 20*log10(abs(PxxA)), 'EdgeColor', 'none');
    axis xy;
    view(0, 90);
    xlim([min(t) max(t)]*1e3);
    ylim([min(f) max(f)]/1e6 + str2double(frequency));
    xlabel('Time [ms]');
    ylabel('Frequency [MHz]');
    ax = gca;
    set(ax, 'clim', [mindB maxdB], 'colormap', turbo);
    c = colorbar;
    c.Label.String = '[dB]';
    c.TickLabelInterpreter = 'latex';
    
%     figure;
%     surf(t*1e3, f/1e6 + str2double(frequency), 20*log10(abs(PxxARFI)), 'EdgeColor', 'none');
%     axis xy;
%     view(0, 90);
%     xlim([min(t) max(t)]*1e3);
%     ylim([min(f) max(f)]/1e6 + str2double(frequency));
%     xlabel('Time [ms]');
%     ylabel('Frequency [MHz]');
%     ax = gca;
%     set(ax, 'clim', [mindB maxdB], 'colormap', turbo);
%     c = colorbar;
%     c.Label.String = '[dB]';
%     c.TickLabelInterpreter = 'latex';
    
    formatFig(gcf, [dataPath 'reconst_results3_' num2str(loopIndex)], 'en', figProp);
    formatFig(gcf, [dataPath 'reconst_results3_' num2str(loopIndex)], 'en', figProp);
%     figure;
%     surf(t*1e3, f/1e6 + str2double(frequency), 20*log10(abs(PxxB)), 'EdgeColor', 'none');
%     axis xy;
%     view(0, 90);
%     xlim([min(t) max(t)]*1e3);
%     ylim([min(f) max(f)]/1e6 + str2double(frequency));
%     xlabel('Time [ms]');
%     ylabel('Frequency [MHz]');
%     ax = gca;
%     set(ax, 'clim', [mindB maxdB], 'colormap', turbo);
%     c = colorbar;
%     c.Label.String = '[dB]';
%     c.TickLabelInterpreter = 'latex';
    close all;  
end

rmpath(['..' filesep '..' filesep 'Misc'])