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
silenceSamples = round(20e-6*params.fs);
totalSamples = numberOfRawSamples;

WinLBlock = [3 19];
SNR = -25;
random_state = 42;

bandwidthVector = 10.72e6;
periodVector = 8.62e-6;

GPSSignals = GPSGen(paramsSignal);
GPSSignals = GPSSignals(1:numberOfRawSamples,:);
GPSSignalsPower = pow_eval(GPSSignals);

monteCarloLoops = 1000;
PfaVector = logspace(-5, 0, 17);
detection_res_cell = cell(length(WinLBlock), 1);

for windowLengthIndex = 1:length(WinLBlock)
    h = window('rectwin', WinLBlock(windowLengthIndex));
    MBlock = fix(totalSamples./WinLBlock(windowLengthIndex));
    detection_res = zeros(monteCarloLoops, MBlock, length(PfaVector));
    
    for Emuindex = 1:monteCarloLoops
        Emuindex
        noise = randn(totalSamples, 1) + 1j*randn(totalSamples, 1);
        noisePower = pow_eval(noise);
        GPSSignalsAux = GPSSignals;
        GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
        mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
        detection_res(Emuindex, :, :) = DeteBlockGoF_FBS(mixtureGPS, h, MBlock, PfaVector);
    end
    
    detection_res_cell{windowLengthIndex} = detection_res;
end

save(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'pai_results' filesep 'results01.mat'], 'detection_res_cell', '-v7.3');

warning('on','all')

rmpath(['..' filesep '..' filesep '.' filesep 'Sigtools' filesep])
rmpath(['..' filesep '..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep '..' filesep  'signalsGeneration' filesep]);
rmpath(['..' filesep '..' filesep  'signalsGeneration' filesep 'sim_params']);