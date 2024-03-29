clear;
clc;
close all;

addpath(['..' filesep '..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep '..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['..' filesep '..' filesep  'signalsGeneration' filesep]);
addpath(['..' filesep '..' filesep  'signalsGeneration' filesep 'sim_params']);

load sim_params_1.mat;

random_state = 42;
rng(random_state)

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
params.JNRVector = 0;
SNR = -25;

numberOfRawSamples = 4096;
totalSamples = numberOfRawSamples;
thresholdVector = -0.3:0.05:0.9;
window_median_length_vector = 0;
monteCarloLoops = 1000;

GPSSignals = GPSGen(paramsSignal);
GPSSignals = GPSSignals(1:numberOfRawSamples,:);
GPSSignalsPower = pow_eval(GPSSignals);
    
detection_res = zeros(monteCarloLoops, length(thresholdVector), length(window_median_length_vector));

for loopIndex = 1:monteCarloLoops
    loopIndex
    noise = randn(totalSamples, 1) + 1j*randn(totalSamples, 1);
    noisePower = pow_eval(noise);
   
    GPSSignalsAux = GPSSignals;
    GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
    mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
    mixtureSignal = mixtureGPS ;
    [W, ~, ~, PxxAux, ~, ~] = nmf_eval_v2(mixtureSignal, params);
%     %------------------------NMF training--------------------------
%     mixtureSignal1 = mixtureSignal(1:round(length(mixtureSignal)/2),:);
%     mixtureSignal2 = mixtureSignal(round(length(mixtureSignal)/2)+1:end,:);
%     [W, ~, ~, ~, ~, ~] = nmf_eval_v2(mixtureSignal1, params);
%     [~, ~, ~, PxxAux, ~, ~] = nmf_eval_v2(mixtureSignal2, params);
%     %--------------------------------------------------------------
%     PxxAux{1, 1} = PxxAux{1, 1}.*W{1, 1};
    inputNMF = abs(PxxAux{1, 1}).^2;
    
    inputNMF = inputNMF - mean(inputNMF);
    inputNMF = inputNMF.*sqrt(1./var(inputNMF));
    inputNMFAux = sqrt(sum(inputNMF.*inputNMF)) + eps;
    inputNMFNormalised = inputNMF./inputNMFAux;
    
    WNormalised = W{1, 1}(:,1) - mean(W{1, 1}(:,1));
%     WNormalised =  W{1, 1}(:,1);
    WNormalised = WNormalised.*sqrt(1./var(WNormalised));
    WNormalised = WNormalised ./ (norm(WNormalised) + eps);
    
    output = inputNMFNormalised.'*WNormalised;
    
    for thresholdIndex = 1:length(thresholdVector)
        for window_median_length_index = 1:length(window_median_length_vector)
            detection_res(loopIndex, thresholdIndex, window_median_length_index) = median(detection_eval(output, thresholdVector(thresholdIndex)));
        end
    end
end

save(['..' filesep '..' filesep '.' filesep 'data' filesep 'TAES_data' filesep 'my_results' filesep 'results44.mat'], 'detection_res', '-v7.3');

rmpath(['..' filesep '..' filesep '.' filesep 'Sigtools' filesep])
rmpath(['..' filesep '..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep '..' filesep  'signalsGeneration' filesep]);
rmpath(['..' filesep '..' filesep  'signalsGeneration' filesep 'sim_params']);