clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['.' filesep 'Sigtools' filesep])
addpath(['.' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['.' filesep 'Misc'])

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;

figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);  

dataPath = ['.' filesep 'figs' filesep '03-27' filesep];

fs = 32.768e6;
numberOfSources = 2;
secondsOfData = 8.62e-6;
secondsOfSilence = 100e-6;
numberOfSamples = secondsOfData*fs;
totalSamples = 4096;
bandwidth = 1e6;
f0 = 0;
random_state = 42;

params.fs = fs;
params.nfft = 1024;
params.nperseg = 128;
params.overlap = params.nperseg/8;
params.numberOfSources = numberOfSources;
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 10000;
params.tolChange = 1e-6;
params.tolError = 1e-6;
params.repetitions = 1;
params.JNRVector = inf;

rng(random_state);

save_fig = false;

%signal mixture definition---------
t = 0:1/fs:(totalSamples/fs - 1/fs);
f = 1e6;
f1 = 1;
f2 = 4;

signal1 = sin(2*pi*f1*f.*t).';
signal2 = sin(2*pi*f2*f.*t).';
signal2Length = length(signal2);

signal1 = [zeros(round(secondsOfSilence*fs), 1); signal1];
signal1Length = length(signal1);

signal2 = [signal2;zeros(signal1Length - signal2Length, 1)];
mixtureSignal = signal1 + signal2;

%--------------------------------------------

[W, H, reconstructError, PxxAux, f, t] = nmf_eval(mixtureSignal, params);

window_length = floor(length(t)/10);

for JNRIndex = 1:length(params.JNRVector)

    for dataIndex = 1:size(PxxAux, 1)

        figure;
        surf(t*1e6, f/1e6, 10*log10(abs(PxxAux{dataIndex, JNRIndex}).^2), 'EdgeColor', 'none');
        axis xy;
        view(0, 90);
        xlim([min(t) max(t)]*1e6);
        ylim([min(f) max(f)]/1e6);
        ylabel('Frequency [MHz]');
        xlabel('Time [$\mu$s]');
        ax = gca;
        set(ax, 'CLim', [-60 20], 'colormap', jet);
        c = colorbar;
        c.Label.String = '[dB]';
        if save_fig
            formatFig(gcf, [dataPath 'spectrogram_sine'], 'en', figProp);
        end


        for i = 1:numberOfSources
            figure;
            plot(f/1e6, W{dataIndex, JNRIndex}(:,i)./ max(W{dataIndex, JNRIndex}(:,i)));
            ylabel('Normalized Magnitude');
            xlabel('Frequency [MHz]');
            ylim([0 1.1])
            xlim([min(f) max(f)]/1e6);

            if save_fig
                formatFig(gcf, [dataPath 'basis_sine_' num2str(i)], 'en', figProp);
            end

            figure;
            plot(t*1e6,  H{dataIndex, JNRIndex}(i,:) ./ max(H{dataIndex, JNRIndex}(i,:)));
            ylabel('Normalized Magnitude');
            xlabel('Time [$\mu$s]');
            ylim([0 1.1])
            xlim([min(t) max(t)]*1e6);
            if save_fig
                formatFig(gcf, [dataPath 'activation_sine_' num2str(i)], 'en', figProp);
            end

            %Mean window
            figure;
            HMean = window_eval(H{dataIndex, JNRIndex}(i,:), window_length, @mean);
            taux = linspace(0, max(t), length(HMean))*1e6;
            plot(taux,  HMean ./ max(HMean));
            ylabel('Normalized Magnitude');
            xlabel('Time [$\mu$s]');
            ylim([0 1.1])
            xlim([min(taux) max(taux)]);
            if save_fig
                formatFig(gcf, [dataPath 'activation_mean_sine_' num2str(i)], 'en', figProp);
            end

        end

        figure;
        idx = find(reconstructError{dataIndex, JNRIndex}, 1, 'last');
        if strcmp(params.betaDivergence, 'frobenius')
            plot(10*log10(reconstructError{dataIndex, JNRIndex}(1:idx)));
        else
            plot(20*log10(reconstructError{dataIndex, JNRIndex}(1:idx)));
        end
        ylabel('Cost Function Error [dB]');
        xlabel('Iterations [$k$]');
        xlim([1 idx])
        if save_fig
            formatFig(gcf, [dataPath 'loss_error_sine'], 'en', figProp);
        end

        for k = 1:numberOfSources
            separatedSignal = W{dataIndex, JNRIndex}(:,k)*H{dataIndex, JNRIndex}(k,:);
            figure;
            surf(t*1e6, f/1e6, 10*log10(separatedSignal + eps(separatedSignal)), 'EdgeColor', 'none');
            ax = gca;
            axis xy;
            view(0, 90);
            xlim([min(t) max(t)]*1e6);
            ylim([min(f) max(f)]/1e6);

            set(ax, 'CLim', [-60 20], 'colormap', jet);
            ylabel('Frequency [MHz]');
            xlabel('Time [$\mu$s]');    
            c = colorbar;
            c.Label.String = '[dB]';
            if save_fig
                formatFig(gcf, [dataPath 'output_spectrogram_sine_' num2str(k)], 'en', figProp);
            end

        end
    end
end

rmpath(['.' filesep 'Misc'])
rmpath(['.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['.' filesep 'Sigtools' filesep])
