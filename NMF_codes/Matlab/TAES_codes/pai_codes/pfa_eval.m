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
totalSamples = numberOfRawSamples;

global Emuindex;                                   % emulate index
global Loopnumb;                                   % loop number
Loopnumb = 1000;

global Segnumb;                                    % number of TF observation intervals within an integration time
Segnumb = 1;
global Nonesegment;                                % number of samples within a TF observation interval
Nonesegment = totalSamples;

global WinLBlock;                                  % window length used in block-wise STFT
WinLBlock = [3 19];
global WinLCano;
WinLCano = 5;
global WintypeCano;
WintypeCano = 1;
global MBlock;                                     % number of samples for each block-wise STFT
MBlock = fix(Nonesegment./WinLBlock);
global Pfa;                                        % false alarm probability for interference detection
Pfa = 1e-4;                                   % false alarm probability for interference detection
global PfaVector
PfaVector = logspace(-8, 0, 17);
global GoFBlockDeteflag;                           % detection flag for GoF-based interference detection algorithm using block-wise STFT
global GoFCanoDeteflag;
JNRVector = 0;

SNR = -25;
random_state = 42;
initialFrequency = 2e6;

bandwidthVector = 10.72e6;
periodVector = 8.62e-6;

global Signal;                                     % aggregated signal
global bandwidthIndex;
global periodIndex;
global windowLengthIndex;

% Run simulations
global GoFBlockDeteflagCell;
GoFBlockDeteflagCell = cell(length(MBlock), 1);

counterLoop = length(bandwidthVector) * length(periodVector) * length(WinLBlock) * Loopnumb;
for bandwidthIndex = 1:length(bandwidthVector)
    for periodIndex = 1:length(periodVector)
        for windowLengthIndex = 1:length(WinLBlock)
            GoFBlockDeteflag = zeros(length(bandwidthVector), length(periodVector), length(PfaVector), MBlock(windowLengthIndex),Segnumb*Loopnumb);

            paramsSignal.Noneperiod = round(periodVector(periodIndex)*params.fs);                   % number of samples with a sweep time
            paramsSignal.IFmin = initialFrequency;                                                  % start frequency
            paramsSignal.IFmax = bandwidthVector(bandwidthIndex) + initialFrequency;                    % end frequency
            paramsSignal.foneperiod(1:paramsSignal.Noneperiod) = linspace(paramsSignal.IFmin, paramsSignal.IFmax, paramsSignal.Noneperiod);
            paramsSignal.Initphase = 0;
            
            [~, GPSSignals] = signalGen(paramsSignal);
            GPSSignals = GPSSignals(1:numberOfRawSamples,:);
            
            GPSSignalsPower = pow_eval(GPSSignals);
            for Emuindex = 1:Loopnumb
                Emuindex
                noise = randn(totalSamples, 1) + 1j*randn(totalSamples, 1);
                noisePower = pow_eval(noise);
                GPSSignalsAux = GPSSignals;
                GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
                mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
                Signal = mixtureGPS;
                DeteBlockGoF
                %         DeteCanoGoF;
            end
            GoFBlockDeteflagCell{windowLengthIndex} = GoFBlockDeteflag;
        end
    end
end

save(['.' filesep 'data' filesep 'resultsPai07.mat'], 'GoFBlockDeteflagCell'); 

warning('on','all')

rmpath(['..' filesep '..' filesep 'aux_files_pai' filesep 'mfiles']);
rmpath(['..' filesep '..' filesep '.' filesep 'Sigtools' filesep])
rmpath(['..' filesep '..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep '..' filesep 'signalsGeneration' filesep 'sim_params']);