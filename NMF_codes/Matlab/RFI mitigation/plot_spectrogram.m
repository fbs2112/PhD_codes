clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '.' filesep 'Misc'])
addpath(['..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['..' filesep  'signalsGeneration' filesep]);

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
dataPath = ['..' filesep '..' filesep '.' filesep 'figs' filesep '2021-05-20' filesep];

load(['..' filesep  'signalsGeneration' filesep 'sim_params' filesep 'sim_params_3.mat']);

params.fs = paramsSignal.Freqsamp;
params.nfft = 256;
params.nperseg = 256;
params.overlap = params.nperseg - 1;
params.hop_size = params.nperseg - params.overlap;
params.window = ones(params.nperseg, 1);

params.window = kaiser(params.nperseg, 60);
% params.window = hamming(params.nperseg);
% params.window = chebwin(params.nperseg, 100);

params.specType = 'power';
params.type = 'power';
params.numberOfSources = 1;
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 10000;
params.tolChange = 1e-6;
params.tolError = 1e-6;
params.repetitions = 1;
SNR = -25;
params.JNRVector = 10;

bandwidthVector = 8e6;
periodVector = (8.62e-6);

initialFrequency = 2e6;

numberOfRawSamples = round(params.fs*1000e-6);

totalSamples = numberOfRawSamples;
monteCarloLoops = 1;

GPSSignals = GPSGen(paramsSignal);
GPSSignals = GPSSignals(1:totalSamples,:);
GPSSignalsPower = pow_eval(GPSSignals);

for loopIndex = 1:monteCarloLoops
    loopIndex
    mixtureSignal = zeros(totalSamples, length(params.JNRVector));
    noise = randn(totalSamples, 1) + 1j*randn(totalSamples, 1);
    noisePower = pow_eval(noise);
    
    for bandwidthIndex = 1:length(bandwidthVector)
        bandwidthIndex
        for periodIndex = 1:length(periodVector)
            paramsSignal.Noneperiod = round(periodVector(periodIndex)*params.fs);                   % number of samples with a sweep time
            paramsSignal.IFmin = initialFrequency;                                                  % start frequency
            paramsSignal.IFmax = bandwidthVector(bandwidthIndex) + initialFrequency;                    % end frequency
            paramsSignal.foneperiod(1:paramsSignal.Noneperiod) = linspace(paramsSignal.IFmin, paramsSignal.IFmax, paramsSignal.Noneperiod);
            paramsSignal.Initphase = 0;
            
            interferenceSignal = interferenceGen(paramsSignal);
            interferenceSignal = interferenceSignal(1:numberOfRawSamples);
            interferenceSignalPower = pow_eval(interferenceSignal);
            
            for i = 1:length(params.JNRVector)
                GPSSignalsAux = GPSSignals;
                interferenceSignalAux = interferenceSignal;
                GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
                mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
                interferenceSignalAux = interferenceSignalAux*sqrt(noisePower*10^(params.JNRVector(i)/10)/interferenceSignalPower);
                mixtureSignal(:,i) = mixtureGPS + interferenceSignalAux;
            end
            
            [W, H, ~, PxxAux, f, t, ~] = nmf_eval_v2(mixtureSignal, params);
            figure;
            surf(t*1e6, f/1e6, 10*log10(abs(PxxAux{1,1}).^2), 'EdgeColor', 'none');
            axis xy;
            view(0, 90);
            xlabel('Time [$\mu$s]');
            ylabel('Frequency [MHz]');
            axis tight;
%             xlim([0 max(t)]*1e6);
%             xticks(0:25:150);
            ylim([min(f) max(f)]/1e6);
            ax = gca;
            set(ax, 'CLim', [0 60], 'colormap', turbo);
            c = colorbar;
            c.Label.String = '[dB]';
            c.TickLabelInterpreter = 'latex';
            
            [Pxx2, f, t] = fsst(mixtureSignal, params.fs, params.window);
            
            figure;
%             mesh(t,f,10*log10(abs(Pxx2)).^2)
%             view(2)
%             formatFig(gcf, [dataPath  'spec_cw_NMF_' num2str(params.nfft) '_' num2str(params.JNRVector)], 'en', figProp);
%             figure;
            surf(t*1e6, f/1e6, 10*log10(abs(Pxx2 + 1e-6).^2), 'EdgeColor', 'none');
            axis xy;
            view(0, 90);
            xlabel('Time [$\mu$s]');
            ylabel('Frequency [MHz]');
            axis tight;
% %             xlim([0 max(t)]*1e6);
% %             xticks(0:25:150);
%             ylim([min(f) max(f)]/1e6);
            ax = gca;
            set(ax, 'CLim', [0 60], 'colormap', turbo);
            c = colorbar;
            c.Label.String = '[dB]';
            c.TickLabelInterpreter = 'latex';    
            
        end
    end
end
rmpath(['..' filesep '.' filesep 'Misc'])
rmpath(['..' filesep '.' filesep 'Sigtools' filesep])
rmpath(['..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep  'signalsGeneration' filesep]);
