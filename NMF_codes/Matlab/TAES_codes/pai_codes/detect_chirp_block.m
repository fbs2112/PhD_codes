clear;
clc;
close all;

warning('off','all')

addpath(['..' filesep '..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep '..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['..' filesep '..' filesep  'signalsGeneration' filesep]);
addpath(['..' filesep '..' filesep  'signalsGeneration' filesep 'sim_params']);

load sim_params_1.mat;

params.fs = paramsSignal.Freqsamp;
numberOfRawSamples = 4096;
totalSamples = numberOfRawSamples;

WinLBlock = 3;
JNRVector = -25:0;
SNR = -25;

initialFrequency = 2e6;
bandwidthVector = (2e6:3e6:14e6);
periodVector = (8.62e-6:1.48e-6:18.97e-6);

GPSSignals = GPSGen(paramsSignal);
GPSSignals = GPSSignals(1:numberOfRawSamples,:);
GPSSignalsPower = pow_eval(GPSSignals);

monteCarloLoops = 100;
PfaVector = logspace(-12, -2, 41);
h = window('rectwin', WinLBlock);
MBlock = fix(totalSamples./WinLBlock);

detection_res = zeros(length(JNRVector), monteCarloLoops, MBlock, length(PfaVector));
pvalue = zeros(length(JNRVector), monteCarloLoops, MBlock);

for Emuindex = 1:monteCarloLoops
    Emuindex
    
    noise = randn(totalSamples, 1) + 1j*randn(totalSamples, 1);
    noisePower = pow_eval(noise);
    GPSSignalsAux = GPSSignals;
    GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
    mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
    
    for JNRIndex = 1:length(JNRVector)
        JNRIndex
        for bandwidthIndex = 1:length(bandwidthVector)
            for periodIndex = 1:length(periodVector)
                paramsSignal.Noneperiod = round(periodVector(periodIndex)*params.fs);                   % number of samples with a sweep time
                paramsSignal.IFmin = initialFrequency;                                     % start frequency
                paramsSignal.IFmax = bandwidthVector(bandwidthIndex) + initialFrequency;                   % end frequency
                paramsSignal.foneperiod(1:paramsSignal.Noneperiod) = linspace(paramsSignal.IFmin, paramsSignal.IFmax, paramsSignal.Noneperiod);
                paramsSignal.Initphase = 0;
                
                interferenceSignal = interferenceGen(paramsSignal);
                interferenceSignal = interferenceSignal(1:numberOfRawSamples);
                interferenceSignalPower = pow_eval(interferenceSignal);
            
                interferenceSignalAux = interferenceSignal;
                interferenceSignalAux = interferenceSignalAux*sqrt(noisePower*10^(JNRVector(JNRIndex)/10)/interferenceSignalPower);
                mixtureSignal = mixtureGPS + interferenceSignalAux;
                [pvalue(JNRIndex, Emuindex, :), detection_res(JNRIndex, Emuindex, :, :)] = DeteBlockGoF_FBS(mixtureSignal, h, MBlock, PfaVector);
            end
        end
    end
end

if isunix
    save(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep 'results_det_9.mat'], 'detection_res', 'pvalue', '-v7.3');
else
    save(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep 'results_det_9.mat'], 'detection_res', 'pvalue', '-v7.3');
end
warning('on','all')

rmpath(['..' filesep '..' filesep '.' filesep 'Sigtools' filesep])
rmpath(['..' filesep '..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep '..' filesep  'signalsGeneration' filesep]);
rmpath(['..' filesep '..' filesep  'signalsGeneration' filesep 'sim_params']);