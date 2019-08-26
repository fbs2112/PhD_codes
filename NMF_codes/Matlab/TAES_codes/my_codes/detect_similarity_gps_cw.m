clear;
clc;
close all;

addpath(['..' filesep '..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep '..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['..' filesep '..' filesep  'signalsGeneration' filesep]);
addpath(['..' filesep '..' filesep  'signalsGeneration' filesep 'sim_params']);

load sim_params_1.mat;

random_state = 42;

params.fs = paramsSignal.Freqsamp;
params.nfft = 512;
params.nperseg = 512;
params.overlap = params.nperseg-1;
params.hop_size = params.nperseg - params.overlap;
params.numberOfSources = 1;
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 10000;
params.tolChange = 1e-6;
params.tolError = 1e-6;
params.repetitions = 1;
SNR = -25;
params.JNRVector = -20:0;

bandwidthVector = 0;
periodVector = 8.72e-6;

rng(random_state)

initialFrequency = params.fs*0.12;
numberOfRawSamples = 4096;
totalSamples = numberOfRawSamples;
thresholdVector = 0:0.005:2;
window_median_length_vector = 51:50:401;
monteCarloLoops = 100;

outputLength = (totalSamples - params.nperseg + 1)/(params.nperseg - params.overlap);
detection_res = zeros(monteCarloLoops, length(bandwidthVector), length(periodVector), ...
    length(params.JNRVector), length(thresholdVector), length(window_median_length_vector), outputLength);

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
            
            GPSSignals = GPSGen(paramsSignal);
            interferenceSignal = interferenceGen(paramsSignal);
            GPSSignals = GPSSignals(1:numberOfRawSamples,:);
            interferenceSignal = interferenceSignal(1:numberOfRawSamples);
            
            interferenceSignalPower = pow_eval(interferenceSignal);
            GPSSignalsPower = pow_eval(GPSSignals);
            
            for i = 1:length(params.JNRVector)
                GPSSignalsAux = GPSSignals;
                interferenceSignalAux = interferenceSignal;
                GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
                mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
                interferenceSignalAux = interferenceSignalAux*sqrt(noisePower*10^(params.JNRVector(i)/10)/interferenceSignalPower);
                mixtureSignal(:,i) = mixtureGPS + interferenceSignalAux;
            end
            
            [W, ~, ~, PxxAux, ~, ~] = nmf_eval_v2(mixtureSignal, params);
            
            for JNRIndex = 1:length(params.JNRVector)
                
                inputNMF = abs(PxxAux{1, JNRIndex}).^2;
                
                inputNMF = inputNMF - mean(inputNMF);
                inputNMF = inputNMF.*sqrt(1./var(inputNMF));
                inputNMFAux = sqrt(sum(inputNMF.*inputNMF)) + eps;
                inputNMFNormalised = inputNMF./inputNMFAux;
                
                WNormalised = W{1, JNRIndex}(:,1) - mean(W{1, JNRIndex}(:,1));
                WNormalised = WNormalised.*sqrt(1./var(WNormalised));
                WNormalised = WNormalised ./ (norm(WNormalised) + eps);
                
                output = inputNMFNormalised.'*WNormalised;
                for thresholdIndex = 1:length(thresholdVector)
                    for window_median_length_index = 1:length(window_median_length_vector)
                        detection_res(loopIndex, bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, window_median_length_index, :) = ...
                            detection_eval(output, thresholdVector(thresholdIndex), window_median_length_vector(window_median_length_index));
                    end
                end
            end
        end
    end
end

if isunix
    save(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' filesep 'results07.mat'], 'detection_res', '-v7.3');
else
end
rmpath(['..' filesep '..' filesep '.' filesep 'Sigtools' filesep])
rmpath(['..' filesep '..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep '..' filesep  'signalsGeneration' filesep]);
rmpath(['..' filesep '..' filesep  'signalsGeneration' filesep 'sim_params']);