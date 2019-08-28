clear;
clc;
close all;

addpath(['..' filesep '..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep '..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['..' filesep '..' filesep  'signalsGeneration' filesep]);
addpath(['..' filesep '..' filesep  'signalsGeneration' filesep 'sim_params']);

load sim_params_1.mat;

params.fs = paramsSignal.Freqsamp;
params.nfft = 64;
params.nperseg = 64;
params.overlap = params.nperseg-1;
params.hop_size = params.nperseg - params.overlap;
params.numberOfSources = 1;
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 10000;
params.tolChange = 1e-6;
params.tolError = 1e-6;
params.repetitions = 1;
params.JNRVector = 0;
SNR = -25;

numberOfRawSamples = 4096;
totalSamples = numberOfRawSamples;
thresholdVector = 0.3:0.005:0.5;
window_median_length_vector = 0;
monteCarloLoops = 1001;

detection_res = zeros(monteCarloLoops, length(thresholdVector), length(window_median_length_vector));

for loopIndex = 1:monteCarloLoops
    loopIndex
    noise = randn(totalSamples, 1) + 1j*randn(totalSamples, 1);
    noisePower = pow_eval(noise);
    
    GPSSignals = GPSGen(paramsSignal);
    GPSSignals = GPSSignals(1:numberOfRawSamples,:);
    GPSSignalsPower = pow_eval(GPSSignals);
    
    GPSSignalsAux = GPSSignals;
    GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
    mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
    mixtureSignal = mixtureGPS ;
    
    if loopIndex == 1
        [W, ~, ~, PxxAux, ~, ~] = nmf_eval_v2(mixtureSignal, params);
        WNormalised = W{1, 1}(:,1) ./ sum(W{1, 1}(:,1));
    else
        [~, ~, ~, PxxAux, ~, ~] = nmf_eval_v2(mixtureSignal, params);
    end
    inputNMF = abs(PxxAux{1, 1}).^2;
    inputNMFNormalised = inputNMF./ sum(inputNMF, 1);
    
    output = sum(inputNMFNormalised .* log(inputNMFNormalised./WNormalised));
    for thresholdIndex = 1:length(thresholdVector)
        for window_median_length_index = 1:length(window_median_length_vector)
            detection_res(loopIndex, thresholdIndex, window_median_length_index) = median(detection_eval(output, thresholdVector(thresholdIndex)));
        end
    end
end

if isunix
    save(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep 'pfa_results' filesep 'results9.mat'], 'detection_res', '-v7.3');
else
    save(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep 'pfa_results' filesep 'results9.mat'], 'detection_res', '-v7.3');
end
rmpath(['..' filesep '..' filesep '.' filesep 'Sigtools' filesep])
rmpath(['..' filesep '..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep '..' filesep  'signalsGeneration' filesep]);
rmpath(['..' filesep '..' filesep  'signalsGeneration' filesep 'sim_params']);