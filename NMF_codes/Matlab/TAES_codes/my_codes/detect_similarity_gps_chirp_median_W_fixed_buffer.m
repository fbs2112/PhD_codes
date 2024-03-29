clear;
clc;
close all;

addpath(['..' filesep '..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep '..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['..' filesep '..' filesep  'signalsGeneration' filesep]);
addpath(['..' filesep '..' filesep  'signalsGeneration' filesep 'sim_params']);

load sim_params_1.mat;
load WChirp.mat;
% 
% WNoise = WNormalised;
random_state = 42;

params.fs = paramsSignal.Freqsamp;
params.nfft = 64;
params.nperseg = 64;
params.overlap = params.nperseg - 1;
params.hop_size = params.nperseg - params.overlap;
params.numberOfSources = 1;
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 10000;
params.tolChange = 1e-6;
params.tolError = 1e-6;
params.repetitions = 1;
SNR = -25;
params.JNRVector = -25:0;
params.JNRVector = 0;

bandwidthVector = 10.72e6;
periodVector = 8.62e-6;

rng(random_state)

initialFrequency = 2e6;
numberOfRawSamples = 4096;
totalSamples = numberOfRawSamples;
thresholdVector = -0.3:0.05:0.9;
window_median_length_vector = 0;
monteCarloLoops = 100;

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
            [~, ~, ~, PxxAux, ~, ~] = nmf_eval_v2(mixtureSignal, params);
%             %------------------------NMF training--------------------------
%             mixtureSignal1 = mixtureSignal(1:round(length(mixtureSignal)/2),:);
%             mixtureSignal2 = mixtureSignal(round(length(mixtureSignal)/2)+1:end,:);
%             [W, ~, ~, ~, ~, ~] = nmf_eval_v2(mixtureSignal1, params);
%             [~, ~, ~, PxxAux, ~, ~] = nmf_eval_v2(mixtureSignal2, params);
%             %--------------------------------------------------------------
            
            for JNRIndex = 1:length(params.JNRVector)
%                 PxxAux{1, JNRIndex} = PxxAux{1, JNRIndex}.*W{1, JNRIndex};
                inputNMF = abs(PxxAux{1, JNRIndex}).^2;
%                 inputNMF = inputNMF.*W{1, JNRIndex};
                
                inputNMF = inputNMF - mean(inputNMF);
                inputNMF = inputNMF.*sqrt(1./var(inputNMF));
                inputNMFAux = sqrt(sum(inputNMF.*inputNMF)) + eps;
                inputNMFNormalised = inputNMF./inputNMFAux;
                
%                 WNormalised = W{1, JNRIndex}(:,1) - mean(W{1, JNRIndex}(:,1));
%                 WNormalised = WNormalised.*sqrt(1./var(WNormalised));
% %                 WNormalised = W{1, JNRIndex};
%                 WNormalised = WNormalised ./ (norm(WNormalised) + eps);
                
                output = inputNMFNormalised.'*WNormalised;
                output = window_eval(output, 100, @mean).';
                for thresholdIndex = 1:length(thresholdVector)
                    for window_median_length_index = 1:length(window_median_length_vector)
                        detection_res(loopIndex, bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, window_median_length_index) = ...
                            median(detection_eval(output, thresholdVector(thresholdIndex)));
                    end
                end
            end
        end
    end
end

save(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results' filesep 'results45.mat'], 'detection_res', '-v7.3');

rmpath(['..' filesep '..' filesep '.' filesep 'Sigtools' filesep])
rmpath(['..' filesep '..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep '..' filesep  'signalsGeneration' filesep]);
rmpath(['..' filesep '..' filesep  'signalsGeneration' filesep 'sim_params']);