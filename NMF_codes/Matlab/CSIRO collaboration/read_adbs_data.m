clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');

addpath(['..' filesep 'Sigtools' filesep])
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);

addpath(['.' filesep 'data']);
load adsb_signal;
load sim_params_2.mat;

nfft = 256;
overlap = nfft - 1;
window = hann(nfft);
fs = 2.4e6;
specType = 'power';
paramsSignal.Freqsamp = fs;
bandwidthVector = 0;
periodVector = 8.62e-6;
initialFrequency = 1e6;

paramsSignal.Noneperiod = round(periodVector*fs);                   % number of samples with a sweep time
paramsSignal.Noneperiod = 51840;                   % number of samples with a sweep time

paramsSignal.IFmin = initialFrequency;                                                  % start frequency
paramsSignal.IFmax = bandwidthVector + initialFrequency;                    % end frequency
paramsSignal.foneperiod(1:paramsSignal.Noneperiod) = linspace(paramsSignal.IFmin, paramsSignal.IFmax, paramsSignal.Noneperiod);
paramsSignal.Initphase = 0;

signalOfInterest = interferenceGen(paramsSignal);

% t = 0:1/fs:16384/fs - 1/fs;
% signalOfInterest = exp(1j*2*pi*initialFrequency.*t);

[PxxAux, f, t] = spectrogram(signalOfInterest, window, overlap, nfft, fs, 'centered', specType);

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

rmpath(['..' filesep 'Sigtools' filesep]);
rmpath(['.' filesep 'data']);
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);