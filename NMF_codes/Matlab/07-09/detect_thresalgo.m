clear;
clc;
close all;

addpath(['..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['..' filesep '.' filesep 'Misc'])

fs = 32.768e6;
numberOfSources = 1;
secondsOfData = 8.62e-6;
secondsOfSilence = 100e-6;
numberOfSamples = secondsOfData*fs;
totalSamples = 4096;
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
params.JNRVector = [-10 -5 0];

rng(random_state);

%signal mixture definition---------
t = 0:1/fs:(secondsOfData - 1/fs);
f = ((bandwidth/2)/secondsOfData)*t + f0;
f1 = 1;

signal1 = exp(1j*2*pi*f1*f.*t).';
signal1 = repmat(signal1, ceil(100e-6/secondsOfData), 1);
onsetTime = 20e-6;
offsetTime = onsetTime + length(signal1)/fs;

signal1 = [zeros(round(onsetTime*fs), 1); signal1;zeros(round(onsetTime*fs), 1)];
signal1Length = length(signal1);

onset = find(signal1, 1, 'first') - params.nperseg;
offset = find(signal1, 1, 'last') - params.nperseg;

mixtureSignal = signal1;
%--------------------------------------------

window_length = round(3e-6*fs);
window_median_length = 51;
similarityName = 'inner';
stdVector = 0;

if strcmp(similarityName, 'gaussian')
    stdVector = 7:2:13;
end

lengthVector = 200:100:500;
stdVector2 = [2 3 4];
forgettingFactorVector = [0.3 0.4 0.5 0.6];
monteCarloLoops = 100;
outputLength = (signal1Length - params.nperseg + 1)/(params.nperseg - params.overlap);

%Pre allocation------------------------------------------------------------
detection_res = zeros(monteCarloLoops, length(params.JNRVector), ...
    length(stdVector), length(lengthVector), length(stdVector2), length(forgettingFactorVector), outputLength);
%--------------------------------------------------------------------------

for loopIndex = 1:monteCarloLoops
    loopIndex
    [W, H, reconstructError, PxxAux, f, t] = nmf_eval(mixtureSignal, params);
    
    for JNRIndex = 1:length(params.JNRVector)
        
        inputNMF = abs(PxxAux{1, JNRIndex}).^2;
        output = zeros(length(t), 1);
        
        for stdIndex = 1:length(stdVector)
            
            for tIndex = 1:length(t)
                if strcmp(similarityName, 'inner')
                    output(tIndex) = similarity_eval(inputNMF(:,tIndex), W{1, JNRIndex}(:,1), similarityName, 'normalized', true);
                else
                    output(tIndex) = similarity_eval(inputNMF(:,tIndex), W{1, JNRIndex}(:,1), similarityName, stdVector(stdIndex), true);
                end
            end
            
            output = TK_filtering(output);
            for lengthVectorIndex = 1:length(lengthVector)
                for stdVector2Index = 1:length(stdVector2)
                    for forgettingFactorVectorIndex = 1:length(forgettingFactorVector)
                        [detection_res(loopIndex, JNRIndex, stdIndex, lengthVectorIndex, ...
                            stdVector2Index, forgettingFactorVectorIndex, :), ~, ~] = ThresholdingAlgo(output, ...
                            lengthVector(lengthVectorIndex), stdVector2(stdVector2Index), forgettingFactorVector(forgettingFactorVectorIndex));
                    end
                end
            end
            
        end
    end
end

save(['..' filesep '.' filesep 'data' filesep '07-09' filesep 'results01.mat'], 'detection_res');

rmpath(['..' filesep '.' filesep 'Misc'])
rmpath(['..' filesep '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep '.' filesep 'Sigtools' filesep])