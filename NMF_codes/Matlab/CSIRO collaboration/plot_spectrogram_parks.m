%This script plots spectrograms to further labelling process
%
%This code was created by Felipe Barboza da Silva
% Copyright (c) School of Engineering. Macquarie University - Australia.
% All rights reserved. 2020.
% DISCLAIMER:
%     This code is strictly private, confidential and personal to its recipients
%     and should not be copied, distributed or reproduced in whole or in part,
%     nor passed to any third party.
% ** Do not duplicate or distribute without written permission from the owners. **

clear;
clc;
close all;

addpath(['..' filesep 'Misc'])

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');

nfft = 1024;
nperseg = 1024;
overlap = nperseg/2;
% window = hann(nperseg);
window = ones(nperseg, 1);

linewidth = 1.5;
fontname = 'Times';
fontsize = 26;
dataPath = [ '..' filesep 'figs' filesep '2021-04-15' filesep];
figProp = struct('size' , fontsize , 'font' , fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);

specType = 'power';
frequency = '1152';
fileName = '2018-09-01-09_52_55_0000000000000000.000000';
diskLetter = 'E:';
year = '2018';

matObj = matfile(['.' filesep 'Labels' filesep year filesep frequency filesep fileName '_labels.mat']);

filePath2 = ['E:' filesep 'Parkes data' filesep year filesep frequency filesep 'signal frames' filesep];

fileIdx = find(matObj.interferenceDetFlag == 1 | matObj.interferenceDetFlag == 3 | matObj.interferenceDetFlag == 5 | matObj.interferenceDetFlag == 6);
% fileIdx = 1:100;
for i = 1:20%length(fileIdx)
    load([filePath2 fileName '_' num2str(fileIdx(i)) '.mat']);
    fs = Fs;
    IF = 24e6;
    t = 0:length(parkesSignalA) - 1;
    carrierConj = exp(-1j*2*pi*(IF)*t/fs).';
%     parkesSignalA = parkesSignalA.*carrierConj;
    figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    
    [Pxx, f, t] = spectrogram(parkesSignalB, window, overlap, nfft, fs, 'centered');
    absPxx = abs(Pxx);
    maxdB = mag2db(max(absPxx(:)));
    
    surf(t*1e3, f/1e6 + str2double(frequency), 20*log10(abs(Pxx)), 'EdgeColor', 'none');
    axis xy;
    view(0, 90);
    xlim([min(t) max(t)]*1e3);
    ylim([min(f)/1e6 max(f)/1e6] + str2double(frequency));
    xlabel('Time [ms]');
    ylabel('Frequency [MHz]');
    ax = gca;
    set(ax, 'colormap', turbo);
    set(ax, 'clim', [20 70], 'colormap', turbo);
    c = colorbar;
    c.Label.String = '[dB]';
    c.TickLabelInterpreter = 'latex';
    close all;
end

rmpath(['..' filesep 'Misc'])