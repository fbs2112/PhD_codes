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

WinLBlock = 19;
JNRVector = 0;
SNR = -25;
random_state = 42;

initialFrequency = params.fs*0.12;
bandwidthVector = 0;
periodVector = 8.62e-6;

paramsSignal.Noneperiod = round(periodVector*params.fs);                   % number of samples with a sweep time
paramsSignal.IFmin = initialFrequency;                                     % start frequency
paramsSignal.IFmax = bandwidthVector + initialFrequency;                   % end frequency
paramsSignal.foneperiod(1:paramsSignal.Noneperiod) = linspace(paramsSignal.IFmin, paramsSignal.IFmax, paramsSignal.Noneperiod);
paramsSignal.Initphase = 0;

GPSSignals = GPSGen(paramsSignal);
GPSSignals = GPSSignals(1:numberOfRawSamples,:);
interferenceSignal = interferenceGen(paramsSignal);
interferenceSignal = interferenceSignal(1:numberOfRawSamples);

GPSSignalsPower = pow_eval(GPSSignals);
interferenceSignalPower = pow_eval(interferenceSignal);

monteCarloLoops = 100;
PfaVector = logspace(-5, 0, 17);
h = window('rectwin', WinLBlock);
MBlock = fix(totalSamples./WinLBlock);

detection_res = zeros(length(JNRVector), monteCarloLoops, MBlock, length(PfaVector));
pvalue = zeros(length(JNRVector), monteCarloLoops, MBlock);

for JNRIndex = 1:length(JNRVector)
   
    for Emuindex = 1:monteCarloLoops
        Emuindex
        noise = randn(totalSamples, 1) + 1j*randn(totalSamples, 1);
        noisePower = pow_eval(noise);
        GPSSignalsAux = GPSSignals;
        interferenceSignalAux = interferenceSignal;
        GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
        mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
        interferenceSignalAux = interferenceSignalAux*sqrt(noisePower*10^(JNRVector(JNRIndex)/10)/interferenceSignalPower);
        mixtureSignal = mixtureGPS + interferenceSignalAux;
        [pvalue(JNRIndex, Emuindex, :), detection_res(JNRIndex, Emuindex, :, :)] = DeteBlockGoF_FBS(mixtureSignal, h, MBlock, PfaVector);
    end
    
end

save(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'pai_results' filesep 'results05.mat'], 'detection_res', '-v7.3');

warning('on','all')

rmpath(['..' filesep '..' filesep '.' filesep 'Sigtools' filesep])
rmpath(['..' filesep '..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep '..' filesep  'signalsGeneration' filesep]);
rmpath(['..' filesep '..' filesep  'signalsGeneration' filesep 'sim_params']);