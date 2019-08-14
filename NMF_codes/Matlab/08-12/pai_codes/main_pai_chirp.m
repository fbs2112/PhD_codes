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
totalSamples = numberOfRawSamples + silenceSamples*2;


global Emuindex;                                   % emulate index
global Loopnumb;                                   % loop number
Loopnumb = 100;

global Segnumb;                                    % number of TF observation intervals within an integration time
Segnumb = 1;
global Nonesegment;                                % number of samples within a TF observation interval
Nonesegment = totalSamples;

global WinLBlock;                                  % window length used in block-wise STFT
WinLBlock = 3;
global WinLCano;
WinLCano = 5;
global WintypeCano;
WintypeCano = 1;
global MBlock;                                     % number of samples for each block-wise STFT
MBlock = fix(Nonesegment/WinLBlock);
global Pfa;                                        % false alarm probability for interference detection
Pfa = 1e-4;                                   % false alarm probability for interference detection
global PfaVector
PfaVector = logspace(-5, 0, 17);
global GoFBlockDeteflag;                           % detection flag for GoF-based interference detection algorithm using block-wise STFT
global GoFCanoDeteflag;
JNRVector = -20:0;
SNR = -25;
random_state = 42;
initialFrequency = 2e6;

bandwidthVector = 10.72e6;
periodVector = 8.62e-6;

global Signal;                                     % aggregated signal
global JNRIndex;
global bandwidthIndex;
global periodIndex;
global counter;
counter = 0;
global counterLoop;
counterLoop = Loopnumb * length(JNRVector) * length(periodVector) * length(bandwidthVector);

% Run simulations

GoFCanoDeteflag = zeros(length(PfaVector), MBlock,Segnumb*Loopnumb);
GoFBlockDeteflag = zeros(length(bandwidthVector), length(periodVector), length(JNRVector), length(PfaVector), MBlock,Segnumb*Loopnumb);

for bandwidthIndex = 1:length(bandwidthVector)
    for periodIndex = 1:length(periodVector)
        
        paramsSignal.Noneperiod = round(periodVector(periodIndex)*params.fs);                   % number of samples with a sweep time
        paramsSignal.IFmin = initialFrequency;                                                  % start frequency
        paramsSignal.IFmax = bandwidthVector(bandwidthIndex) + initialFrequency;                    % end frequency
        paramsSignal.foneperiod(1:paramsSignal.Noneperiod) = linspace(paramsSignal.IFmin, paramsSignal.IFmax, paramsSignal.Noneperiod);
        paramsSignal.Initphase = 0;

        [interferenceSignal, GPSSignals] = signalGen(paramsSignal);
        GPSSignals = GPSSignals(1:numberOfRawSamples,:);
        interferenceSignal = interferenceSignal(1:numberOfRawSamples);

        GPSSignals = [zeros(silenceSamples, size(GPSSignals, 2)); GPSSignals; zeros(silenceSamples, size(GPSSignals, 2))];
        interferenceSignal = [zeros(silenceSamples, 1); interferenceSignal; zeros(silenceSamples, 1)];

        interferenceSignalPower = pow_eval(interferenceSignal);
        GPSSignalsPower = pow_eval(GPSSignals);


        for JNRIndex = 1:length(JNRVector)
            JNRIndex
            
            for Emuindex = 1:Loopnumb
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
                counter = counter + 1;
        %         DeteCanoGoF;
            end
        end
    end
end
warning('on','all')

rmpath(['..' filesep '..' filesep 'aux_files_pai' filesep 'mfiles']);
rmpath(['..' filesep '..' filesep '.' filesep 'Sigtools' filesep])
rmpath(['..' filesep '..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep '..' filesep 'signalsGeneration' filesep 'sim_params']);