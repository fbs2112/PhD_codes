%This script trains the NMF according to a given ADS-B signal and the
%signal of interest (CW)
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
    
addpath(['..' filesep 'Sigtools' filesep])
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms']);
addpath(['.' filesep 'data']);

load ADSB_label.mat;

trueLabels = find(interferenceDetFlag);
falseLabels = find(~interferenceDetFlag);

numberOfSignalFrames = 10;

interferenceFrames = zeros(1280000, numberOfSignalFrames);
nonInterferenceFrames = zeros(1280000, numberOfSignalFrames);

monteCarloLoops = 1;
params.JNRVector = 0;
params.fs = 128e6;
params.nfft = 256;
params.nperseg = 256;
params.overlap = params.nperseg/2;
params.window = ones(params.nperseg, 1);
params.hop_size = params.nperseg - params.overlap;
params.specType = 'power';
params.numberOfSources = 10;
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 500;
params.tolChange = 1e-6;
params.tolError = 1e-6;
params.repetitions = 1;
params.verbose = 1;
JNR = 20;

freqVector = linspace(-64, -56, params.nfft)*1e6;

for i = 1:numberOfSignalFrames
    load(['.' filesep 'data' filesep 'dataParks_' num2str(trueLabels(i)) '.mat']);
    interferenceFrames(:,i) = parksSignal;
    
%     [Pxx, f, t] = spectrogram(parksSignal, params.window, params.overlap, params.nfft, params.fs, params.specType);
%     figure;
%     surf(t*1e3, f/1e6, 10*log10(abs(Pxx)), 'EdgeColor', 'none');
%     axis xy;
%     view(0, 90);
%     xlim([min(t) max(t)]*1e3);
% %     ylim([min(f)/1e6 -56]);
%     xlabel('Time [ms]');
%     ylabel('Frequency [MHz]');
%     ax = gca;
%     set(ax, 'colormap', jet);
%     c = colorbar;
%     c.Label.String = '[dB]';
%     c.TickLabelInterpreter = 'latex';
    
    load(['.' filesep 'data' filesep 'dataParks_' num2str(falseLabels(i)) '.mat']);
    nonInterferenceFrames(:,i) = parksSignal;
    
%     [Pxx, f, t] = spectrogram(parksSignalAux, params.window, params.overlap, freqVector, params.fs, params.specType);
%     figure;
%     surf(t*1e3, f/1e6, 10*log10(abs(Pxx)), 'EdgeColor', 'none');
%     axis xy;
%     view(0, 90);
%     xlim([min(t) max(t)]*1e3);
%     ylim([min(f)/1e6 -56]);
%     xlabel('Time [ms]');
%     ylabel('Frequency [MHz]');
%     ax = gca;
%     set(ax, 'colormap', jet);
%     c = colorbar;
%     c.Label.String = '[dB]';
%     c.TickLabelInterpreter = 'latex';
end

interferenceFrames = interferenceFrames(:);
nonInterferenceFrames = nonInterferenceFrames(:);

for loopIndex = 1:monteCarloLoops
   
    [WInterf, HInterf, errorInterfTrain, PxxInterf, ~, ~] = nmf_eval_v2(interferenceFrames, params);                %ADS-B frequency content evaluation
    [WSignal, HSignal, errorSignalTrain, PxxSignal, ~, ~] = nmf_eval_v2(nonInterferenceFrames, params);  %CW frequency content evaluation
    
    W0 = [WInterf{1,1} WSignal{1,1}];
end

save(['.' filesep 'data' filesep 'nmf_training_ADSB_01.mat'], 'W0');             %saving trained matrix W

rmpath(['.' filesep 'data']);
rmpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms']);
rmpath(['..' filesep 'Sigtools' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep]);