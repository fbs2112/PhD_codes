clear;
clc;
close all;

warning('off','all')

addpath(['..' filesep '..' filesep 'aux_files_pai' filesep 'mfiles']);
addpath(['..' filesep '..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep '..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep '..' filesep 'signalsGeneration' filesep 'sim_params']);

load sim_params_1.mat;

params.fs = paramsSignal.Freqsamp;
numberOfRawSamples = 4096;
silenceSamples = round(20e-6*params.fs);
% totalSamples = numberOfRawSamples + silenceSamples*2;
totalSamples = numberOfRawSamples;


global Emuindex;                                   % emulate index
global Loopnumb;                                   % loop number
Loopnumb = 1000;

global Segnumb;                                    % number of TF observation intervals within an integration time
Segnumb = 1;
global Nonesegment;                                % number of samples within a TF observation interval
Nonesegment = totalSamples;

global WinLBlock;                                  % window length used in block-wise STFT
WinLBlock = 19;
global WinLCano;
WinLCano = 5;
global WintypeCano;
WintypeCano = 1;
global MBlock;                                     % number of samples for each block-wise STFT
MBlock = fix(Nonesegment/WinLBlock);
global Pfa;                                        % false alarm probability for interference detection
Pfa = 1e-4;                                   % false alarm probability for interference detection
global PfaVector
PfaVector = logspace(-5, 0, 25);
global GoFBlockDeteflag;                           % detection flag for GoF-based interference detection algorithm using block-wise STFT
global GoFCanoDeteflag;
% JNRVector = -20:0;
JNRVector = -17;


SNR = -25;
random_state = 42;
initialFrequency = params.fs*0.12;

bandwidthVector = 0;
periodVector = 8.72e-6;

paramsSignal.Noneperiod = round(periodVector*params.fs);                   % number of samples with a sweep time
paramsSignal.IFmin = initialFrequency;                                                  % start frequency
paramsSignal.IFmax = bandwidthVector + initialFrequency;                    % end frequency
paramsSignal.foneperiod(1:paramsSignal.Noneperiod) = linspace(paramsSignal.IFmin, paramsSignal.IFmax, paramsSignal.Noneperiod);
paramsSignal.Initphase = 0;

[interferenceSignal, GPSSignals] = signalGen(paramsSignal);
GPSSignals = GPSSignals(1:numberOfRawSamples,:);
interferenceSignal = interferenceSignal(1:numberOfRawSamples);

% GPSSignals = [zeros(silenceSamples, size(GPSSignals, 2)); GPSSignals; zeros(silenceSamples, size(GPSSignals, 2))];
% interferenceSignal = [zeros(silenceSamples, 1); interferenceSignal; zeros(silenceSamples, 1)];

interferenceSignalPower = pow_eval(interferenceSignal);
GPSSignalsPower = pow_eval(GPSSignals);

global Signal;                                     % aggregated signal
global JNRIndex;
% Run simulations

for JNRIndex = 1:length(JNRVector)
    JNRIndex
    GoFBlockDeteflag = zeros(length(PfaVector), MBlock,Segnumb*Loopnumb);
    GoFCanoDeteflag = zeros(length(PfaVector), MBlock,Segnumb*Loopnumb);
    for Emuindex = 1:Loopnumb
        Emuindex
        noise = randn(totalSamples, 1) + 1j*randn(totalSamples, 1);
        noisePower = pow_eval(noise);
        GPSSignalsAux = GPSSignals;
        interferenceSignalAux = interferenceSignal;
        GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
        mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
        interferenceSignalAux = interferenceSignalAux*sqrt(noisePower*10^(JNRVector(JNRIndex)/10)/interferenceSignalPower);
        mixtureSignal = mixtureGPS + interferenceSignalAux;
        Signal = mixtureSignal;

        DeteBlockGoF;       
        DeteCanoGoF;
    end
end
warning('on','all')

rmpath(['..' filesep '..' filesep 'aux_files_pai' filesep 'mfiles']);
rmpath(['..' filesep '..' filesep '.' filesep 'Sigtools' filesep])
rmpath(['..' filesep '..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep '..' filesep 'signalsGeneration' filesep 'sim_params']);