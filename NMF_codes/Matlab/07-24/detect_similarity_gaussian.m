clear;
clc;
close all;

addpath(['..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep '.' filesep 'Sigtools' filesep 'NMF_algorithms'])

fs = 32.768e6;
numberOfSources = 1;
secondsOfData = 8.62e-6;
bandwidth = 1e6;
f0 = 0;
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
params.JNRVector = [-20 -15 -10 -5 0];

rng(random_state);

%signal mixture definition---------
t = 0:1/fs:(secondsOfData - 1/fs);
f = ((bandwidth/2)/secondsOfData)*t + f0;
f1 = 1;

signal1 = exp(1j*2*pi*f1*f.*t).';
signal1 = repmat(signal1, ceil(100e-6/secondsOfData), 1);
onsetTime = 20e-6;

signal1 = [zeros(round(onsetTime*fs), 1); signal1;zeros(round(onsetTime*fs), 1)];
mixtureSignal = signal1;
%--------------------------------------------

similarityName = 'gaussian';
stdVector = 0;

if strcmp(similarityName, 'gaussian')
    stdVector = 7:13;
end

thresholdVector = 0.1:0.05:0.9;
window_median_length_vector = 51:50:401;
monteCarloLoops = 100;

outputLength = (length(signal1) - params.nperseg + 1)/(params.nperseg - params.overlap);
detection_res = zeros(monteCarloLoops, length(params.JNRVector), length(thresholdVector), length(stdVector), length(window_median_length_vector), outputLength);

for loopIndex = 1:monteCarloLoops
    loopIndex
    [W, ~, ~, PxxAux, ~, ~] = nmf_eval(mixtureSignal, params);
    
    for JNRIndex = 1:length(params.JNRVector)
        
        for thresholdIndex = 1:length(thresholdVector)
            
            inputNMF = abs(PxxAux{1, JNRIndex}).^2;
            for stdIndex = 1:length(stdVector)
                output = zeros(outputLength, 1);
                for tIndex = 1:outputLength
                    if strcmp(similarityName, 'inner')
                        output(tIndex) = similarity_eval(inputNMF(:,tIndex), W{1, JNRIndex}(:,1), similarityName, 'normalized', true);
                    elseif strcmp(similarityName, 'gaussian')
                        output(tIndex) = similarity_eval(inputNMF(:,tIndex), W{1, JNRIndex}(:,1), similarityName, stdVector(stdIndex), true);
                    else
                        output(tIndex) = similarity_eval(inputNMF(:,tIndex), W{1, JNRIndex}(:,1), similarityName, [], false);
                    end
                end
                
                outputNormalised = output ./ max(output);
                
                for window_median_length_index = 1:length(window_median_length_vector)
                    detection_res(loopIndex, JNRIndex, thresholdIndex, stdIndex, window_median_length_index, :) = detection_eval(outputNormalised, thresholdVector(thresholdIndex), ...
                        window_median_length_vector(window_median_length_index));
                end
            end
        end
    end
end

save(['..' filesep '.' filesep 'data' filesep '07-24' filesep 'results04.mat'], 'detection_res', '-v7.3');

rmpath(['..' filesep '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep '.' filesep 'Sigtools' filesep])
exit;