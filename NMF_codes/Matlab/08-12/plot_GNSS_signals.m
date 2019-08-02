clear;
clc;
% close all

set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); 
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaulttextInterpreter','latex')


addpath(['..' filesep '.' filesep 'data' filesep '08-12']);  

load gnss_signals.mat

fs = 32.768e6;
params.fs = fs;
params.nfft = 128;
params.nperseg = 128;
params.overlap = params.nperseg-1;
params.hop_size = params.nperseg - params.overlap;

Signal = Signal(1:4096);

if isreal(Signal)
    [PxxAux, f, t] = spectrogram(Signal, hann(params.nperseg), params.overlap, params.nfft, params.fs, 'power');
else
    [PxxAux, f, t] = spectrogram(Signal, hann(params.nperseg), params.overlap, params.nfft, params.fs, 'centered', 'power');
end


figure;
surf(t*1e6, f/1e6, 10*log10(abs(PxxAux).^2), 'EdgeColor', 'none');
axis xy;
view(0, 90);
xlim([min(t) max(t)]*1e6);
ylim([min(f) max(f)]/1e6);
ylabel('Frequency [MHz]');
xlabel('Time [$\mu$s]');
ax = gca;
set(ax, 'CLim', [-10 50], 'colormap', jet);
c = colorbar;
c.Label.String = '[dB]';

Signalass = Signalass(1:4096);

if isreal(Signal)
    [PxxAux, f, t] = spectrogram(Signalass, hann(params.nperseg), params.overlap, params.nfft, params.fs, 'power');
else
    [PxxAux, f, t] = spectrogram(Signalass, hann(params.nperseg), params.overlap, params.nfft, params.fs, 'centered', 'power');
end

figure;
surf(t*1e6, f/1e6, 10*log10(abs(PxxAux).^2), 'EdgeColor', 'none');
axis xy;
view(0, 90);
xlim([min(t) max(t)]*1e6);
ylim([min(f) max(f)]/1e6);
ylabel('Frequency [MHz]');
xlabel('Time [$\mu$s]');
ax = gca;
set(ax, 'CLim', [-10 50], 'colormap', jet);
c = colorbar;
c.Label.String = '[dB]';

aux1 = Signal;
aux2 = Signalass;

load gnss_signals2.mat;


Signal = Signal(1:4096);

if isreal(Signal)
    [PxxAux, f, t] = spectrogram(Signal, hann(params.nperseg), params.overlap, params.nfft, params.fs, 'power');
else
    [PxxAux, f, t] = spectrogram(Signal, hann(params.nperseg), params.overlap, params.nfft, params.fs, 'centered', 'power');
end


figure;
surf(t*1e6, f/1e6, 10*log10(abs(PxxAux).^2), 'EdgeColor', 'none');
axis xy;
view(0, 90);
xlim([min(t) max(t)]*1e6);
ylim([min(f) max(f)]/1e6);
ylabel('Frequency [MHz]');
xlabel('Time [$\mu$s]');
ax = gca;
set(ax, 'CLim', [-10 50], 'colormap', jet);
c = colorbar;
c.Label.String = '[dB]';

Signalass = Signalass(1:4096);

if isreal(Signal)
    [PxxAux, f, t] = spectrogram(Signalass, hann(params.nperseg), params.overlap, params.nfft, params.fs, 'power');
else
    [PxxAux, f, t] = spectrogram(Signalass, hann(params.nperseg), params.overlap, params.nfft, params.fs, 'centered', 'power');
end

figure;
surf(t*1e6, f/1e6, 10*log10(abs(PxxAux).^2), 'EdgeColor', 'none');
axis xy;
view(0, 90);
xlim([min(t) max(t)]*1e6);
ylim([min(f) max(f)]/1e6);
ylabel('Frequency [MHz]');
xlabel('Time [$\mu$s]');
ax = gca;
set(ax, 'CLim', [-10 50], 'colormap', jet);
c = colorbar;
c.Label.String = '[dB]';

rmpath(['..' filesep '.' filesep 'data' filesep '08-12']);  

addpath(['..' filesep 'Sigtools' filesep])
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params' filesep]);

load sim_params_1.mat;

SNR = -25;
JNR = 0;




