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

numberOfTrainingCells = 20;
numberOfGuardCells = 4;
rng(random_state);

%signal mixture definition---------
t = 0:1/fs:(secondsOfData - 1/fs);
f = ((bandwidth/2)/secondsOfData)*t + f0;
f1 = 1;
numberOfTrials = 100;

signal1 = exp(1j*2*pi*f1);
noise = 1/sqrt(2)*(randn(length(t), numberOfTrials) + 1j*randn(length(t), numberOfTrials));

idxs = randperm(length(t), 4);

noise(idxs,:) = signal1;

data = abs(noise).^2;
%--------------------------------------------

 detector = phased.CFARDetector('NumTrainingCells', numberOfTrainingCells, 'NumGuardCells', numberOfGuardCells,...
     'ProbabilityFalseAlarm', 0.01, 'ThresholdOutputPort', true);
 %             detector.ThresholdFactor = 'Auto';
 detector.Method = 'CA';
 detector.Rank = 5;
 
[detection_res, ~] = detector(data, 1:length(data));

idxsSignal = zeros(length(t), 1);
idxsSignal(idxs) = 1;

pd = mean(detection_res(idxs,:), 2);

aux = 1:length(t);
aux(idxs) = 0;
aux2 = aux(aux~=0);

pfa = mean(mean(detection_res(aux2,:), 2));

rmpath(['..' filesep '.' filesep 'Sigtools' filesep])
rmpath(['..' filesep '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
%%

clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep '.' filesep 'Sigtools' filesep 'NMF_algorithms'])

fs = 32.768e6;
bandwidth = 1e6;
f0 = 0;
secondsOfData = 8.62e-6;
random_state = 42;

JNR = 0;
numberOfTrainingCells = 10;
numberOfGuardCells = 4;
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
mixtureSignal = signal1;
powMixture = pow_eval(mixtureSignal);
signal1Length = length(mixtureSignal);

noise = randn(signal1Length, 1);

if ~isreal(mixtureSignal)
    noise = randn(signal1Length, 1) + 1j*randn(signal1Length, 1);
end

powNoise = pow_eval(noise);
powAux = powMixture/db2pow(JNR);
noise = noise*sqrt(powAux/powNoise);
data = (abs(mixtureSignal + noise).^2);

%--------------------------------------------

 detector = phased.CFARDetector('NumTrainingCells', numberOfTrainingCells, 'NumGuardCells', numberOfGuardCells,...
     'ProbabilityFalseAlarm', 0.01, 'ThresholdOutputPort', true);
 %             detector.ThresholdFactor = 'Auto';
 detector.Method = 'CA';
 detector.Rank = 5;
 
[detection_res, thres] = detector(data, 1:length(data));

detection_res_median = window_eval(detection_res, 5, @median);

figure;
plot(data)
hold on;
plot(detection_res);
figure;
plot(detection_res_median);

figure;
plot(thres);
rmpath(['..' filesep '.' filesep 'Sigtools' filesep])
rmpath(['..' filesep '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
