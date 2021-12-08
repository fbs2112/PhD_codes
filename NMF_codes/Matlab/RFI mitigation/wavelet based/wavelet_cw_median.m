clear;
clc;
close all;

addpath(['..' filesep '..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep '..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['..' filesep '..' filesep  'signalsGeneration' filesep]);
addpath(['..' filesep '..' filesep  'signalsGeneration' filesep 'sim_params']);

load sim_params_1.mat;

params.fs = paramsSignal.Freqsamp;
params.nfft = 16;
params.nperseg = 16;
params.overlap = params.nperseg - 1;
params.hop_size = params.nperseg - params.overlap;
params.window = ones(params.nperseg, 1);
params.specType = 'power';
params.type = 'power';
params.numberOfSources = 1;
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 10000;
params.tolChange = 1e-6;
params.tolError = 1e-6;
params.repetitions = 1;
SNR = -25;
params.JNRVector = -25:0;

bandwidthVector = 0;
periodVector = 8.62e-6;

initialFrequency = params.fs*0.12;
numberOfRawSamples = 4096;
totalSamples = numberOfRawSamples;
thresholdVector = 0:0.05:10;
window_median_length_vector = 0;
monteCarloLoops = 1000;

outputLength = (totalSamples - params.nperseg + 1)/(params.nperseg - params.overlap);
detection_res = zeros(monteCarloLoops, length(bandwidthVector), length(periodVector), ...
    length(params.JNRVector), length(thresholdVector), length(window_median_length_vector));

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
            
            for i = 1:length(params.JNRVector)
                GPSSignalsAux = GPSSignals;
                interferenceSignalAux = interferenceSignal;
                GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
                mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
                interferenceSignalAux = interferenceSignalAux*sqrt(noisePower*10^(params.JNRVector(i)/10)/interferenceSignalPower);
                mixtureSignal(:,i) = mixtureGPS + interferenceSignalAux;
            end
                        
           for JNRIndex = 1:length(params.JNRVector)
                wav_output = wav_eval(mixtureSignal(:,JNRIndex));
                
                for thresholdIndex = 1:length(thresholdVector)
                    for window_median_length_index = 1:length(window_median_length_vector)
                        detection_res_aux = zeros(16, 1);
                        for wav_Index = 1:16
                            detection_res_aux(wav_Index) = median(detection_eval(wav_output(:,wav_Index), thresholdVector(thresholdIndex)));
                        end
                        detection_res(loopIndex, bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, window_median_length_index) = any(detection_res_aux);
                    end
                end
            end
        end
    end
end

save(['.' filesep 'data' filesep 'results_det_23.mat'], 'detection_res', '-v7.3');

rmpath(['..' filesep '..' filesep '.' filesep 'Sigtools' filesep])
rmpath(['..' filesep '..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep '..' filesep  'signalsGeneration' filesep]);
rmpath(['..' filesep '..' filesep  'signalsGeneration' filesep 'sim_params']);