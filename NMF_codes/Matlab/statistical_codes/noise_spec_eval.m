clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['..' filesep  'signalsGeneration' filesep]);
addpath(['..' filesep  'signalsGeneration' filesep 'sim_params']);

load sim_params_1.mat;

warning('off','all')

params.fs = paramsSignal.Freqsamp;
params.nfft = 128;
params.nperseg = 128;
params.overlap = params.nperseg - 1;
params.hop_size = params.nperseg - params.overlap;
params.window = ones(params.nperseg, 1);
params.specType = 'power';
params.numberOfSources = 1;
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 10000;
params.tolChange = 1e-6;
params.tolError = 1e-6;
params.repetitions = 1;
JNRVector = -10;
SNR = -25;


bandwidthVector = 2e6;
periodVector = 8.62e-6;
initialFrequency = 2e6;

numberOfRawSamples = 4096;
totalSamples = numberOfRawSamples;
thresholdVector = 0:0.05:0.5;
monteCarloLoops = 100;

timeBins = floor((totalSamples - params.overlap)/(params.nperseg - params.overlap));
detection_res = zeros(monteCarloLoops, timeBins, length(thresholdVector));
pvalue = zeros(monteCarloLoops, timeBins);
kvalue = zeros(monteCarloLoops, timeBins);

% for loopIndex = 1:monteCarloLoops
%     loopIndex
    noise = randn(totalSamples, 1) + 1j*randn(totalSamples, 1);
    noisePower = pow_eval(noise);
    
    GPSSignals = GPSGen(paramsSignal);
    GPSSignals = GPSSignals(1:numberOfRawSamples,:);
    GPSSignalsPower = pow_eval(GPSSignals);
    
    GPSSignalsAux = GPSSignals;
    GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
    mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
    mixtureSignal = mixtureGPS ;
    
    [PxxAux, f, t] = spectrogram(mixtureSignal, params.window, params.overlap, params.nfft, params.fs, 'centered', params.specType);
    spec = abs(PxxAux).^2;
    
    figure;
    surf(t*1e6, f/1e6, 10*log10(spec), 'EdgeColor', 'none');
    axis xy;
    axis tight;
    colormap(jet); 
    view(0,90);
    ylabel('Frequency [MHz]');
    xlabel('Time [$\mu$ s]');
    colorbar;
    caxis([0 40]);
    
    averageSpecNoise = mean(spec, 2);
   
    figure;
    plot(10*log10(fftshift(abs(fft(averageSpecNoise)))))
    
    noise = randn(totalSamples, 1) + 1j*randn(totalSamples, 1);
    noisePower = pow_eval(noise);
    
    paramsSignal.Noneperiod = round(periodVector*params.fs);                   % number of samples with a sweep time
    paramsSignal.IFmin = initialFrequency;                                                  % start frequency
    paramsSignal.IFmax = bandwidthVector + initialFrequency;                    % end frequency
    paramsSignal.foneperiod(1:paramsSignal.Noneperiod) = linspace(paramsSignal.IFmin, paramsSignal.IFmax, paramsSignal.Noneperiod);
    paramsSignal.Initphase = 0;
    
    interferenceSignal = interferenceGen(paramsSignal);
    interferenceSignal = interferenceSignal(1:numberOfRawSamples);
    interferenceSignalPower = pow_eval(interferenceSignal);
   
    GPSSignalsAux = GPSSignals;
    interferenceSignalAux = interferenceSignal;
    GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
    mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
    interferenceSignalAux = interferenceSignalAux*sqrt(noisePower*10^(JNRVector/10)/interferenceSignalPower);
    mixtureSignal = mixtureGPS + interferenceSignalAux;
    
    [PxxAux, ~, ~] = spectrogram(mixtureSignal, params.window, params.overlap, params.nfft, params.fs, 'centered', params.specType);
    specSignal = abs(PxxAux).^2;
    figure;
    surf(t*1e6, f/1e6, 10*log10(specSignal), 'EdgeColor', 'none');
    axis xy;
    axis tight;
    colormap(jet); 
    view(0,90);
    ylabel('Frequency [MHz]');
    xlabel('Time [$\mu$ s]');
    colorbar;
    caxis([0 40]);
    
    specSignalFilt = max(specSignal - averageSpecNoise, 0);
    
    figure;
    surf(t*1e6, f/1e6, 10*log10(specSignalFilt), 'EdgeColor', 'none');
    axis xy;
    axis tight;
    colormap(jet); 
    view(0,90);
    ylabel('Frequency [MHz]');
    xlabel('Time [$\mu$ s]');
    colorbar;
    caxis([0 40]);
% end

% save(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
%         'Doctorate' filesep 'Research' filesep 'data' filesep 'statistical_data' filesep ...
%         filesep 'pfa_results' filesep 'results1.mat'], 'detection_res', '-v7.3');

rmpath(['..' filesep '.' filesep 'Sigtools' filesep])
rmpath(['..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep  'signalsGeneration' filesep]);
rmpath(['..' filesep  'signalsGeneration' filesep 'sim_params']);
