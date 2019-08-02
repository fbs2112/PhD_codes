clear;
clc;
close all;

addpath(['..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);

load sim_params_1.mat;

fs = paramsSignal.Freqsamp;

numberOfSources = 1;
random_state = 42;

params.fs = fs;
params.nfft = 128;
params.nperseg = 128;
params.overlap = params.nperseg-1;
params.hop_size = params.nperseg - params.overlap;
params.numberOfSources = numberOfSources;
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 10000;
params.tolChange = 1e-6;
params.tolError = 1e-6;
params.repetitions = 1;
SNR = -25;
params.JNRVector = [-10];

rng(random_state)

numberOfRawSamples = 4096;
silenceSamples = round(20e-6*fs);

[interferenceSignal, GPSSignals] = signalGen(paramsSignal);
GPSSignals = GPSSignals(1:numberOfRawSamples,:);
interferenceSignal = interferenceSignal(1:numberOfRawSamples);

GPSSignals = [zeros(silenceSamples, size(GPSSignals, 2)); GPSSignals; zeros(silenceSamples, size(GPSSignals, 2))];
interferenceSignal = [zeros(silenceSamples, 1); interferenceSignal; zeros(silenceSamples, 1)];

onset = find(GPSSignals(:,1), 1, 'first') - params.nperseg;
offset = find(GPSSignals(:,1), 1, 'last') - params.nperseg;

interferenceSignalPower = pow_eval(interferenceSignal);
GPSSignalsPower = pow_eval(GPSSignals);

thresholdVector = 0.1:0.05:0.9;
window_median_length_vector = 51:50:401;
monteCarloLoops = 100;

outputLength = (numberOfRawSamples + silenceSamples*2 - params.nperseg + 1)/(params.nperseg - params.overlap);
detection_res = zeros(monteCarloLoops, length(params.JNRVector), length(thresholdVector), length(window_median_length_vector), outputLength);

for loopIndex = 1:monteCarloLoops
    loopIndex
    mixtureSignal = zeros(numberOfRawSamples + silenceSamples*2, length(params.JNRVector));
    noise = randn(numberOfRawSamples + silenceSamples*2, 1) + 1j*randn(numberOfRawSamples + silenceSamples*2, 1);
    noisePower = pow_eval(noise);
    
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
        outputNormalised = output ./ max(output);
        figure;
        plot(output);
       
        for thresholdIndex = 1:length(thresholdVector)
            for window_median_length_index = 1:length(window_median_length_vector)
                detection_res(loopIndex, JNRIndex, thresholdIndex, window_median_length_index, :) = detection_eval(outputNormalised, thresholdVector(thresholdIndex), ...
                    window_median_length_vector(window_median_length_index));
            end
        end
    end
end

% save(['..' filesep '.' filesep 'data' filesep '08-12' filesep 'results02.mat'], 'detection_res', '-v7.3');

rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep '.' filesep 'Sigtools' filesep])