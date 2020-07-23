%This script plots the NMF output spectrograms as well as the training
%ADS-B signal's spectrogram and the received signal's spectrogram
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

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');

addpath(['.' filesep 'data']);

load adsb_signal.mat;
load nmf_testing_03.mat;

nfft = 256;
overlap = nfft - 1;
window = ones(nfft, 1);
fs = 2.4e6;
specType = 'power';

%Training ADS-B signal's spectrogram
rcvTrain = rcv(1.5e4:1.75e4);
[PxxAux, f, t] = spectrogram(rcvTrain, window, overlap, nfft, fs, 'centered', specType);

figure;
surf(t*1e3, f/1e6, 10*log10(abs(PxxAux)), 'EdgeColor', 'none');
axis xy;
view(0, 90);
xlim([min(t) max(t)]*1e3);
ylim([min(f) max(f)]/1e6);
ylabel('Frequency [MHz]');
xlabel('Time [ms]');
ax = gca;
set(ax, 'CLim', [-5 10], 'colormap', jet);
c = colorbar;
c.Label.String = '[dB]';
c.TickLabelInterpreter = 'latex';

%Testing ADS-B signal's spectrogram
[PxxAux, f, t] = spectrogram(rcv(end/2+1:end), window, overlap, nfft, fs, 'centered', specType);
figure;
surf(t*1e3, f/1e6, 10*log10(abs(PxxAux)), 'EdgeColor', 'none');
axis xy;
view(0, 90);
xlim([min(t) max(t)]*1e3);
ylim([min(f) max(f)]/1e6);
ylabel('Frequency [MHz]');
xlabel('Time [ms]');
ax = gca;
set(ax, 'CLim', [-5 10], 'colormap', jet);
c = colorbar;
c.Label.String = '[dB]';
c.TickLabelInterpreter = 'latex';

%Received mixture signal's spectrogram
[PxxAux, f, t] = spectrogram(mixtureSignal, window, overlap, nfft, fs, 'centered', specType);
figure;
surf(t*1e3, f/1e6, 10*log10(abs(PxxAux)), 'EdgeColor', 'none');
axis xy;
view(0, 90);
xlim([min(t) max(t)]*1e3);
ylim([min(f) max(f)]/1e6);
ylabel('Frequency [MHz]');
xlabel('Time [ms]');
ax = gca;
set(ax, 'CLim', [-15 40], 'colormap', jet);
c = colorbar;
c.Label.String = '[dB]';
c.TickLabelInterpreter = 'latex';
close all;
%ADS-B signal's reconstructed spectrogram
figure;
surf(t*1e3, f/1e6, 10*log10(abs(S(:,:,1))), 'EdgeColor', 'none');
axis xy;
view(0, 90);
xlim([min(t) max(t)]*1e3);
ylim([min(f) max(f)]/1e6);
ylabel('Frequency [MHz]');
xlabel('Time [ms]');
ax = gca;
set(ax, 'CLim', [-15 40], 'colormap', jet);
c = colorbar;
c.Label.String = '[dB]';
c.TickLabelInterpreter = 'latex';

%CW's reconstructed spectrogram
figure;
surf(t*1e3, f/1e6, 10*log10(abs(S(:,:,2))), 'EdgeColor', 'none');
axis xy;
view(0, 90);
xlim([min(t) max(t)]*1e3);
ylim([min(f) max(f)]/1e6);
ylabel('Frequency [MHz]');
xlabel('Time [ms]');
ax = gca;
set(ax, 'CLim', [-15 40], 'colormap', jet);
c = colorbar;
c.Label.String = '[dB]';
c.TickLabelInterpreter = 'latex';

rmpath(['.' filesep 'data']);