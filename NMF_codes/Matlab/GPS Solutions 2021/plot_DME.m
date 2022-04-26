clear;
clc;
close all;

addpath(['..' filesep 'Misc'])

linewidth = 1.5;
fontname = 'Times';
fontsize = 26;
dataPath = ['..' filesep 'figs' filesep '2021-04-15' filesep];
figProp = struct('size' , fontsize , 'font' , fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);

matObj = matfile(['.' filesep 'DME signals' filesep 'signals_12MHz_train_FL390.mat']);

fs = 12e6;
nfft = 1024;
nperseg = 1024;
overlap = nperseg/2;
hop_size = nperseg - overlap;
window = ones(nperseg, 1);

[Pxx, f, t] = spectrogram(matObj.DME(1,1:fs*10e-3,1), window, overlap, nfft, fs, 'centered');

figure;
surf(t*1e3, f/1e6, 20*log10(abs(Pxx) + 1e-8), 'EdgeColor', 'none');
axis xy;
view(0, 90);
xlim([min(t) max(t)]*1e3);
xlabel('Time [ms]');
ylabel('Frequency [MHz]');
ax = gca;
set(ax, 'colormap', turbo);
% set(ax, 'clim', [10 90], 'colormap', turbo);
c = colorbar;
c.Label.String = '[dB]';
c.TickLabelInterpreter = 'latex';

formatFig(gcf, [dataPath 'DME_train'], 'en', figProp);

matObj2 = matfile(['.' filesep 'DME signals' filesep 'signals_12MHz_test_FL390.mat']);

[Pxx, f, t] = spectrogram(matObj2.DME(1,1:fs*10e-3), window, overlap, nfft, fs, 'centered');

figure;
surf(t*1e3, f/1e6, 20*log10(abs(Pxx) + 1e-8), 'EdgeColor', 'none');
axis xy;
view(0, 90);
xlim([min(t) max(t)]*1e3);
xlabel('Time [ms]');
ylabel('Frequency [MHz]');
ax = gca;
set(ax, 'colormap', turbo);
% set(ax, 'clim', [10 90], 'colormap', turbo);
c = colorbar;
c.Label.String = '[dB]';
c.TickLabelInterpreter = 'latex';
formatFig(gcf, [dataPath 'DME_test'], 'en', figProp);


figure;
[h,w] = pspectrum(matObj.DME, fs);
plot(w/1e6 + 1176,10*log10(abs(h)));
xlabel('Frequency [MHz]');
ylabel('Power Spectrum [dB]');
axis tight
formatFig(gcf, [dataPath 'pow_spectrum'], 'en', figProp);

rmpath(['..' filesep 'Misc'])
