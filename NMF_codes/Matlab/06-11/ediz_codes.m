clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['..' filesep '.' filesep 'Misc'])

fs = 32.768e6;
numberOfSources = 1;
secondsOfData = 8.62e-6;
secondsOfSilence = 100e-6;
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
params.JNRVector = [-20 -15 -10 -5 0 10];

params.JNRVector = inf;


rng(random_state);


signalLength = 125e-6;
numberOfSamples = round(signalLength*fs);

%signal mixture definition---------
signal = randn(numberOfSamples, 1) + 1j*randn(numberOfSamples, 1);
signalPower = signal'*signal/numberOfSamples;

desiredSignalPower = db2pow(10);
signal = signal*sqrt(desiredSignalPower/signalPower);
mixtureSignal = signal;

signalLength = length(signal);
%--------------------------------------------

% %signal mixture definition---------
% t = 0:1/fs:(secondsOfData - 1/fs);
% f = ((bandwidth/2)/secondsOfData)*t + f0;
% f1 = 1;
% 
% signal1 = exp(1j*2*pi*f1*f.*t).';
% signal1 = repmat(signal1, ceil(100e-6/secondsOfData), 1);
% 
% onsetTime = 40e-6;
% offsetTime = onsetTime + length(signal1)/fs;
% signal1 = [zeros(round(onsetTime*fs), 1); signal1; zeros(round(40e-6*fs), 1)];
% signal1Length = length(signal1);
% mixtureSignal = signal1;
% %--------------------------------------------

window_length = round(3e-6*fs);
similarityName = 'inner';
stdVector = 0;

if strcmp(similarityName, 'gaussian')
    stdVector = 7:2:13;
end

thresholdVector = 1;
monteCarloLoops = 50;

output = zeros((signalLength - params.nperseg + 1)/(params.nperseg - params.overlap), length(params.JNRVector), monteCarloLoops);
outputVar = zeros((signalLength - params.nperseg + 1)/(params.nperseg - params.overlap), length(params.JNRVector), monteCarloLoops);
outputTK = zeros((signalLength - params.nperseg + 1)/(params.nperseg - params.overlap), length(params.JNRVector), monteCarloLoops);
outputTKVar = zeros((signalLength - params.nperseg + 1)/(params.nperseg - params.overlap), length(params.JNRVector), monteCarloLoops);


for loopIndex = 1:monteCarloLoops
    loopIndex
    [W, H, reconstructError, PxxAux, f, t] = nmf_eval(mixtureSignal, params);
    
    for JNRIndex = 1:length(params.JNRVector)
        
        inputNMF = abs(PxxAux{1, JNRIndex}).^2;
        
        
        for stdIndex = 1:length(stdVector)
            
            for tIndex = 1:length(t)
                if strcmp(similarityName, 'inner')
                    output(tIndex, JNRIndex, loopIndex) = similarity_eval(inputNMF(:,tIndex), W{1, JNRIndex}(:,1), similarityName, 'normalized', true);
                else
                    output(tIndex, JNRIndex, loopIndex) = similarity_eval(inputNMF(:,tIndex), W{1, JNRIndex}(:,1), similarityName, stdVector(stdIndex), true);
                end
            end     
            
            outputVar(:, JNRIndex, loopIndex) = window_eval(output(:, JNRIndex, loopIndex), window_length, @var);
            outputTK(:, JNRIndex, loopIndex) = TK_filtering(output(:, JNRIndex, loopIndex));
            outputTKVar(:, JNRIndex, loopIndex) = window_eval(outputTK(:, JNRIndex, loopIndex), window_length, @var);
           
        end
    end
end

save(['..' filesep '.' filesep 'data' filesep '06-11' filesep 'resultsEdiz2.mat'], 'output', 'outputVar', 'outputTK', 'outputTKVar');

rmpath(['..' filesep '.' filesep 'Misc'])
rmpath(['..' filesep '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep '.' filesep 'Sigtools' filesep])