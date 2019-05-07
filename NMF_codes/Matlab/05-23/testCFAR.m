clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

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

% params.JNRVector = [-20 -15 -10 -5 0];
params.JNRVector = [30];
numberOfTrainingCells = 200;
numberOfGuardCells = 20;
rng(random_state);

%signal mixture definition---------
t = 0:1/fs:(secondsOfData - 1/fs);
f = ((bandwidth/2)/secondsOfData)*t + f0;
f1 = 1;

signal1 = exp(1j*2*pi*f1*f.*t).';
signal1 = repmat(signal1, ceil(100e-6/secondsOfData), 1);
onsetTime = 200e-6;
offsetTime = onsetTime + length(signal1)/fs;

signal1 = [zeros(round(onsetTime*fs), 1); signal1;zeros(round(onsetTime*fs), 1)];

mixtureSignal = signal1;

powMixture = pow_eval(mixtureSignal);
signal1Length = length(mixtureSignal);

noise = randn(signal1Length, 1);

if ~isreal(mixtureSignal)
    noise = randn(signal1Length, 1) + 1j*randn(signal1Length, 1);
end
powNoise = pow_eval(noise);
powAux = powMixture/db2pow(params.JNRVector);
noise2 = noise*sqrt(powAux/powNoise);
data = mixtureSignal + noise2;
data = abs(data).^2;
%--------------------------------------------

 detector = phased.CFARDetector('NumTrainingCells', numberOfTrainingCells, 'NumGuardCells', numberOfGuardCells,...
     'ProbabilityFalseAlarm', 0.01, 'ThresholdOutputPort', true);
 %             detector.ThresholdFactor = 'Auto';
 detector.Method = 'CA';
 detector.Rank = 5;
 %             detector.ThresholdFactor = 'Custom';
 %             detector.CustomThresholdFactor = 0.5;

 
 outputMeanBuffer = buffer(data, numberOfTrainingCells+numberOfGuardCells + 1, numberOfTrainingCells+numberOfGuardCells);
 [detection_res, thres] = detector(outputMeanBuffer, (numberOfTrainingCells+numberOfGuardCells)/2 + 1);
 
%  [detection_res, thres] = detector(data, 1:length(data));
 figure
 plot(data);
 hold on
 plot(detection_res);
 figure;
 plot(thres)

rmpath(['..' filesep '.' filesep 'Sigtools' filesep])
rmpath(['..' filesep '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
