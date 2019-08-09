clear;
clc;
close all;

warning('off','all')

addpath(['..' filesep '..' filesep 'aux_files_pai' filesep 'mfiles']);
addpath(['..' filesep '..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep '..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep '..' filesep 'signalsGeneration' filesep 'sim_params']);

load sim_params_1.mat;

global Emuindex;                                   % emulate index
global Loopnumb;                                   % loop number
Loopnumb = 100;

global Segnumb;                                    % number of TF observation intervals within an integration time
Segnumb = 8;
global Nonesegment;                                % number of samples within a TF observation interval
Nonesegment = 4096;

global WinLBlock;                                  % window length used in block-wise STFT
WinLBlock = 19;
global MBlock;                                     % number of samples for each block-wise STFT
MBlock = fix(Nonesegment/WinLBlock);
global Pfa;                                        % false alarm probability for interference detection
Pfa = 1e-4;                                   % false alarm probability for interference detection
global PfaVector
PfaVector = [1e-5 1e-4 1e-3 1e-2 1e-1 1];
global GoFBlockDeteflag;                           % detection flag for GoF-based interference detection algorithm using block-wise STFT
GoFBlockDeteflag = zeros(length(PfaVector), MBlock,Segnumb*Loopnumb);

params.fs = paramsSignal.Freqsamp;
JNRVector = 20;
SNR = -25;
random_state = 42;
initialFrequency = 2e6;
numberOfRawSamples = 4096;
bandwidthVector = 2e6;
periodVector = 8.72e-6;

paramsSignal.Noneperiod = round(periodVector*params.fs);                   % number of samples with a sweep time
paramsSignal.IFmin = initialFrequency;                                                  % start frequency
paramsSignal.IFmax = bandwidthVector + initialFrequency;                    % end frequency
paramsSignal.foneperiod(1:paramsSignal.Noneperiod) = linspace(paramsSignal.IFmin, paramsSignal.IFmax, paramsSignal.Noneperiod);
paramsSignal.Initphase = 0;

mixtureSignal = zeros(numberOfRawSamples*Segnumb, length(JNRVector));
noise = randn(numberOfRawSamples*Segnumb, 1) + 1j*randn(numberOfRawSamples*Segnumb, 1);
noisePower = pow_eval(noise);

[~, GPSSignals] = signalGen(paramsSignal);
t = 0:1/params.fs:(numberOfRawSamples*Segnumb/params.fs - 1/params.fs);
f = params.fs*0.12;
interferenceSignal = sin(2*pi*f.*t).';
interferenceSignalPower = pow_eval(interferenceSignal);
GPSSignalsPower = pow_eval(GPSSignals);

for i = 1:length(JNRVector)
    GPSSignalsAux = GPSSignals;
    interferenceSignalAux = interferenceSignal;
    GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
    mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
    interferenceSignalAux = interferenceSignalAux*sqrt(noisePower*10^(JNRVector(i)/10)/interferenceSignalPower);
    mixtureSignal(:,i) = mixtureGPS + interferenceSignalAux;
end

global Signal;                                     % aggregated signal
Signal = mixtureSignal;

% Run simulations
for Emuindex = 1:Loopnumb   
    Emuindex
    DeteBlockGoF;   
end
warning('on','all')

rmpath(['..' filesep '..' filesep 'aux_files_pai' filesep 'mfiles']);
rmpath(['..' filesep '..' filesep '.' filesep 'Sigtools' filesep])
rmpath(['..' filesep '..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep '..' filesep 'signalsGeneration' filesep 'sim_params']);