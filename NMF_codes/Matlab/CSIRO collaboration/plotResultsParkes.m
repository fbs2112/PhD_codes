clear;
clc;
close all;

addpath(['.' filesep 'data']);

load ADSB_label.mat;
load nmf_testing_ADSB_01.mat;

trueLabels = find(interferenceDetFlag);
falseLabels = find(~interferenceDetFlag);

numberOfSignalFrames = 10;
monteCarloLoops = numberOfSignalFrames;

params.fs = 128e6;
params.nfft = 256;
params.nperseg = 256;
params.overlap = params.nperseg/2;
params.hop_size = params.nperseg - params.overlap;
params.window = ones(params.nperseg, 1);
params.specType = 'power';

interferenceFrames = zeros(1280000, numberOfSignalFrames);
nonInterferenceFrames = zeros(1280000, numberOfSignalFrames);

for i = 1:numberOfSignalFrames
    load(['.' filesep 'data' filesep 'dataParks_' num2str(trueLabels(i + numberOfSignalFrames)) '.mat']);
    
    [Pxx, f, t] = spectrogram(parksSignal, params.window, params.overlap, params.nfft, params.fs, params.specType);
    figure;
    surf(t*1e3, f/1e6, 10*log10(abs(Pxx)), 'EdgeColor', 'none');
    axis xy;
    view(0, 90);
    xlim([min(t) max(t)]*1e3);
    ylim([min(f) max(f)]/1e6);
    xlabel('Time [ms]');
    ylabel('Frequency [MHz]');
    ax = gca;
    set(ax, 'clim', [-5 40], 'colormap', jet);
    c = colorbar;
    c.Label.String = '[dB]';
    c.TickLabelInterpreter = 'latex';
    
    [Pxx, f, t] = spectrogram(xHat(:,1,1,i), params.window, params.overlap, params.nfft, params.fs, params.specType);
    figure;
    surf(t*1e3, f/1e6, 10*log10(abs(Pxx)), 'EdgeColor', 'none');
    axis xy;
    view(0, 90);
    xlim([min(t) max(t)]*1e3);
    ylim([min(f) max(f)]/1e6);
    xlabel('Time [ms]');
    ylabel('Frequency [MHz]');
    ax = gca;
    set(ax, 'clim', [-5 40], 'colormap', jet);
    c = colorbar;
    c.Label.String = '[dB]';
    c.TickLabelInterpreter = 'latex';
    
    [Pxx, f, t] = spectrogram(xHat(:,2,1,i), params.window, params.overlap, params.nfft, params.fs, params.specType);
    figure;
    surf(t*1e3, f/1e6, 10*log10(abs(Pxx)), 'EdgeColor', 'none');
    axis xy;
    view(0, 90);
    xlim([min(t) max(t)]*1e3);
    ylim([min(f) max(f)]/1e6);
    xlabel('Time [ms]');
    ylabel('Frequency [MHz]');
    ax = gca;
    set(ax, 'clim', [-5 40], 'colormap', jet);
    c = colorbar;
    c.Label.String = '[dB]';
    c.TickLabelInterpreter = 'latex';
    close all;
end

rmpath(['.' filesep 'data']);