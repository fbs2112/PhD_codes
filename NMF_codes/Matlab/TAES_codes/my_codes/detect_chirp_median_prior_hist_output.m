clear;
clc;
close all;

addpath(['..' filesep '..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep '..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['..' filesep '..' filesep  'signalsGeneration' filesep]);
addpath(['..' filesep '..' filesep  'signalsGeneration' filesep 'sim_params']);

load sim_params_1.mat;
load WNoise.mat;

WNoise = WNormalised;
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
% params.JNRVector = 0;

bandwidthVector = 10.72e6;
periodVector = 8.62e-6;

rng(random_state)

initialFrequency = 2e6;
numberOfRawSamples = 4096;
totalSamples = numberOfRawSamples;
thresholdVector = -0.1:0.005:0.2;
window_median_length_vector = 0;
monteCarloLoops = 100;

detection_res = zeros(monteCarloLoops, length(bandwidthVector), length(periodVector), ...
    length(params.JNRVector), length(thresholdVector), length(window_median_length_vector));

GPSSignals = GPSGen(paramsSignal);
GPSSignals = GPSSignals(1:numberOfRawSamples,:);
GPSSignalsPower = pow_eval(GPSSignals);

output_res = zeros(monteCarloLoops, length(params.JNRVector));
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
            if loopIndex == 1
                [W, ~, ~, PxxAux, ~, ~] = nmf_eval_v2(mixtureSignal, params);
            else
                [~, ~, ~, PxxAux, ~, ~] = nmf_eval_v2(mixtureSignal, params);
            end
            
            for JNRIndex = 1:length(params.JNRVector)
                
                inputNMF = abs(PxxAux{1, JNRIndex}).^2;
                
                inputNMF = inputNMF - mean(inputNMF);
                inputNMF = inputNMF.*sqrt(1./var(inputNMF));
                inputNMFAux = sqrt(sum(inputNMF.*inputNMF)) + eps;
                inputNMFNormalised = inputNMF./inputNMFAux;
                
                output = inputNMFNormalised.'*WNoise;
                output_res(loopIndex, JNRIndex) = median(output);
             
            end
        end
    end
end

if isunix
    save(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep 'results_det_3.mat'], 'detection_res', '-v7.3');
else
    save(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep 'results_det_3.mat'], 'detection_res', '-v7.3');
end

rmpath(['..' filesep '..' filesep '.' filesep 'Sigtools' filesep])
rmpath(['..' filesep '..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep '..' filesep  'signalsGeneration' filesep]);
rmpath(['..' filesep '..' filesep  'signalsGeneration' filesep 'sim_params']);