clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');

addpath(['.' filesep 'data']);
load adsb_signal.mat;
load nmf_testing_01.mat;

nfft = 256;
overlap = nfft - 1;
window = ones(nfft, 1);
fs = 2.4e6;
specType = 'power';

[PxxAux, f, t] = spectrogram(rcv(1:end/2), window, overlap, nfft, fs, 'centered', specType);

figure;
surf(t*1e6, f/1e6, 10*log10(abs(PxxAux)), 'EdgeColor', 'none');
axis xy;
view(0, 90);
xlim([min(t) max(t)]*1e6);
ylim([min(f) max(f)]/1e6);
ylabel('Frequency [MHz]');
xlabel('Time [$\mu$s]');
ax = gca;
set(ax, 'CLim', [-5 10], 'colormap', jet);
set(ax, 'colormap', jet);

c = colorbar;
c.Label.String = '[dB]';
c.TickLabelInterpreter = 'latex';

[PxxAux, f, t] = spectrogram(rcv(end/2+1:end), window, overlap, nfft, fs, 'centered', specType);

figure;
surf(t*1e6, f/1e6, 10*log10(abs(PxxAux)), 'EdgeColor', 'none');
axis xy;
view(0, 90);
xlim([min(t) max(t)]*1e6);
ylim([min(f) max(f)]/1e6);
ylabel('Frequency [MHz]');
xlabel('Time [$\mu$s]');
ax = gca;
set(ax, 'CLim', [-5 10], 'colormap', jet);
set(ax, 'colormap', jet);

c = colorbar;
c.Label.String = '[dB]';
c.TickLabelInterpreter = 'latex';

[PxxAux, f, t] = spectrogram(mixtureSignal, window, overlap, nfft, fs, 'centered', specType);
figure;
surf(t*1e6, f/1e6, 10*log10(abs(PxxAux)), 'EdgeColor', 'none');
axis xy;
view(0, 90);
xlim([min(t) max(t)]*1e6);
ylim([min(f) max(f)]/1e6);
ylabel('Frequency [MHz]');
xlabel('Time [$\mu$s]');
ax = gca;
set(ax, 'CLim', [-15 40], 'colormap', jet);
set(ax, 'colormap', jet);
c = colorbar;
c.Label.String = '[dB]';
c.TickLabelInterpreter = 'latex';

figure;
surf(t*1e6, f/1e6, 10*log10(abs(S(:,:,1))), 'EdgeColor', 'none');
axis xy;
view(0, 90);
xlim([min(t) max(t)]*1e6);
ylim([min(f) max(f)]/1e6);
ylabel('Frequency [MHz]');
xlabel('Time [$\mu$s]');
ax = gca;
set(ax, 'CLim', [-15 40], 'colormap', jet);
set(ax, 'colormap', jet);
c = colorbar;
c.Label.String = '[dB]';
c.TickLabelInterpreter = 'latex';

figure;
surf(t*1e6, f/1e6, 10*log10(abs(S(:,:,2))), 'EdgeColor', 'none');
axis xy;
view(0, 90);
xlim([min(t) max(t)]*1e6);
ylim([min(f) max(f)]/1e6);
ylabel('Frequency [MHz]');
xlabel('Time [$\mu$s]');
ax = gca;
set(ax, 'CLim', [-15 40], 'colormap', jet);
set(ax, 'colormap', jet);
c = colorbar;
c.Label.String = '[dB]';
c.TickLabelInterpreter = 'latex';

rmpath(['.' filesep 'data']);