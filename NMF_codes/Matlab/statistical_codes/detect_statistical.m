clear;
clc;
close all;

addpath(['..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['..' filesep  'signalsGeneration' filesep]);
addpath(['..' filesep  'signalsGeneration' filesep 'sim_params']);

load sim_params_1.mat;
warning('off','all')

params.fs = paramsSignal.Freqsamp;
params.nfft = 1024;
params.nperseg = 1024;
params.overlap = params.nperseg/2;
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
SNR = -25;
params.JNRVector = -25:0;

% initialFrequency = params.fs*0.12;
initialFrequency = 2e6;

bandwidthVector = 3e6;
periodVector = 8.62e-6;
numberOfRawSamples = 4096;
totalSamples = numberOfRawSamples;
thresholdVector = 0:0.05:0.5;
monteCarloLoops = 100;

timeBins = floor((totalSamples - params.overlap)/(params.nperseg - params.overlap));
detection_res = zeros(monteCarloLoops, length(bandwidthVector), length(periodVector), ...
    length(params.JNRVector), timeBins, length(thresholdVector));
pvalue = zeros(monteCarloLoops, length(bandwidthVector), length(periodVector), ...
    length(params.JNRVector), timeBins);

GPSSignals = GPSGen(paramsSignal);
GPSSignals = GPSSignals(1:numberOfRawSamples,:);
GPSSignalsPower = pow_eval(GPSSignals);

for loopIndex = 1:monteCarloLoops
    loopIndex
    mixtureSignal = zeros(totalSamples, length(params.JNRVector));
    noise = randn(totalSamples, 1) + 1j*randn(totalSamples, 1);
    noisePower = pow_eval(noise);
    
    for bandwidthIndex = 1:length(bandwidthVector)
        for periodIndex = 1:length(periodVector)
            paramsSignal.Noneperiod = round(periodVector(periodIndex)*params.fs);                   % number of samples with a sweep time
            paramsSignal.IFmin = initialFrequency;                                                  % start frequency
            paramsSignal.IFmax = bandwidthVector(bandwidthIndex) + initialFrequency;                    % end frequency
            paramsSignal.foneperiod(1:paramsSignal.Noneperiod) = linspace(paramsSignal.IFmin, paramsSignal.IFmax, paramsSignal.Noneperiod);
            paramsSignal.Initphase = 0;
           
            interferenceSignal = interferenceGen(paramsSignal);
            interferenceSignal = interferenceSignal(1:numberOfRawSamples);
            interferenceSignalPower = pow_eval(interferenceSignal);
            
            for JNRIndex = 1:length(params.JNRVector)
                
                GPSSignalsAux = GPSSignals;
                interferenceSignalAux = interferenceSignal;
                GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
                mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
                interferenceSignalAux = interferenceSignalAux*sqrt(noisePower*10^(params.JNRVector(JNRIndex)/10)/interferenceSignalPower);
                mixtureSignal = mixtureGPS + interferenceSignalAux;
                
                [detection_res(loopIndex, bandwidthIndex, periodIndex, JNRIndex, :, :), pvalue(loopIndex, bandwidthIndex, periodIndex, JNRIndex, :), ~] = ...
                    deteTimeSlice(mixtureSignal, params, thresholdVector);
             
            end
        end
    end
end

save(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'statistical_data' filesep 'results_det_2.mat'], 'detection_res', '-v7.3');

rmpath(['..' filesep '.' filesep 'Sigtools' filesep])
rmpath(['..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep  'signalsGeneration' filesep]);
rmpath(['..' filesep  'signalsGeneration' filesep 'sim_params']);
warning('on','all')