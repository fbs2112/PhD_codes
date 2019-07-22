clear;
clc;
close all;

addpath(['..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep '.' filesep 'Sigtools' filesep 'NMF_algorithms'])

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
params.JNRVector = [-20 -15 -10 -5 0];

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

similarityName = 'inner';
stdVector = 0;

if strcmp(similarityName, 'gaussian')
    stdVector = 7:2:13;
end

thresholdVector = 1;
window_length_vector = 50:50:200;
window_median_length_vector = 51:50:401;
monteCarloLoops = 100;

numberOfTrainingCellsVector = 200:100:800;
numberOfGuardCellsVector = 10:10:30;

outputLength = (length(signal1) - params.nperseg + 1)/(params.nperseg - params.overlap);
detection_res = zeros(monteCarloLoops, length(params.JNRVector), length(thresholdVector), length(window_length_vector),...
    length(numberOfTrainingCellsVector), length(numberOfGuardCellsVector), length(window_median_length_vector), outputLength);

detector = phased.CFARDetector('NumTrainingCells', numberOfTrainingCellsVector(1), 'NumGuardCells', numberOfGuardCellsVector(1),...
    'ProbabilityFalseAlarm', 0.05, 'ThresholdOutputPort', true);
detector.Method = 'GOCA';

for loopIndex = 1:monteCarloLoops
    loopIndex
    [W, ~, ~, PxxAux, ~, ~] = nmf_eval(mixtureSignal, params);
    
    for JNRIndex = 1:length(params.JNRVector)
        
        for thresholdIndex = 1:length(thresholdVector)
            
            inputNMF = abs(PxxAux{1, JNRIndex}).^2;
            for stdIndex = 1:length(stdVector)
                
                inputNMF = inputNMF - mean(inputNMF);
                inputNMF = inputNMF.*sqrt(1./var(inputNMF));
                inputNMFAux = sqrt(sum(inputNMF.*inputNMF)) + eps;
                inputNMFNormalised = inputNMF./inputNMFAux;
                
                WNormalised = W{1, JNRIndex}(:,1) - mean(W{1, JNRIndex}(:,1));
                WNormalised = WNormalised.*sqrt(1./var(WNormalised));
                WNormalised = WNormalised ./ (norm(WNormalised) + eps);

                output = inputNMFNormalised.'*WNormalised;
                for window_length_index = 1:length(window_length_vector)
                    output = window_eval(output, window_length_vector(window_length_index), @var);
                    
                    for numberOfTrainingCellsIndex = 1:length(numberOfTrainingCellsVector)
                        detector.NumTrainingCells = numberOfTrainingCellsVector(numberOfTrainingCellsIndex);
                        for numberOfGuardCellsIndex = 1:length(numberOfGuardCellsVector)
                            detector.NumGuardCells = numberOfGuardCellsVector(numberOfGuardCellsIndex);
                            for window_median_length_index = 1:length(window_median_length_vector)
                                [detection_res_aux, ~] = detector(output.', 1:length(output));
                                detection_res(loopIndex, JNRIndex, thresholdIndex, window_length_index,...
                                     numberOfTrainingCellsIndex, numberOfGuardCellsIndex, window_median_length_index, :) = window_eval(detection_res_aux, window_median_length_vector(window_median_length_index), @median);
                            end
                            release(detector);
                        end
                    end
                end
            end
        end
    end
end

save(['..' filesep '.' filesep 'data' filesep '07-24' filesep 'results05.mat'], 'detection_res', '-v7.3');

rmpath(['..' filesep '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep '.' filesep 'Sigtools' filesep])