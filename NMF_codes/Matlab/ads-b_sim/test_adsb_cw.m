clear;
clc;
close all;

addpath(['..' filesep 'Sigtools' filesep])
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
addpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['.' filesep 'data']);

load nmf_training_03.mat;
load adsb_signal.mat;
load sim_params_2.mat;

monteCarloLoops = 1;
SNR = -25;
params.JNRVector = [10];

params.fs = 2.4e6;
params.nfft = 256;
params.nperseg = 256;
params.overlap = params.nperseg - 1;
params.hop_size = params.nperseg - params.overlap;
params.window = ones(params.nperseg, 1);
params.specType = 'power';
params.numberOfSources = 10;
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 500;
params.tolChange = 1e-6;
params.tolError = 1e-6;
params.repetitions = 1;

JNRVector = params.JNRVector;

% rcv = rcv_data(:,2);
rcv = rcv(end/2+1:end);
totalSamples = length(rcv);
fs = 2.4e6;
bandwidthVector = 0;
periodVector = 8.62e-6;
initialFrequency = 1e6;

params.numberOfSources = params.numberOfSources*2;
params.init = 'custom';
params.W0 = W0;
params.transform = false;
params.verbose = true;

interferencePower = pow_eval(rcv);

mixtureSignal = single(zeros(totalSamples, length(params.JNRVector), monteCarloLoops));
xHat = zeros(length(mixtureSignal), 2, length(params.JNRVector), monteCarloLoops);
paramsSignal.Intenumb = length(rcv);

for loopIndex = 1:monteCarloLoops
    loopIndex
    paramsSignal.Freqsamp = fs;  
    paramsSignal.Noneperiod = round(periodVector*fs);                   % number of samples with a sweep time
    paramsSignal.IFmin = initialFrequency;                                                  % start frequency
    paramsSignal.IFmax = bandwidthVector + initialFrequency;                    % end frequency
    paramsSignal.foneperiod(1:paramsSignal.Noneperiod) = linspace(paramsSignal.IFmin, paramsSignal.IFmax, paramsSignal.Noneperiod);
    paramsSignal.Initphase = 0;
    
    signalOfInterest = interferenceGen(paramsSignal);
    signalOfInterest = signalOfInterest(1:length(rcv));
    signalOfInterestPower = pow_eval(signalOfInterest);
    
    for JNRIndex = 1:length(params.JNRVector)
        interferenceSignalAux = rcv;
        interferenceMultiplier = sqrt(signalOfInterestPower*10.^(JNRVector(JNRIndex)/10)./interferencePower);
        mixtureSignal(:,JNRIndex,loopIndex) = (interferenceSignalAux.*interferenceMultiplier) + single(signalOfInterest);
    end
    
    [WTest, HTest, errorTest, PxxTest, f, t] = nmf_eval_v2(mixtureSignal(:,:,loopIndex), params);
    S = zeros([size(PxxTest{1,1}) 2]);
    
    for JNRIndex = 1:length(params.JNRVector)
        for i = 1:2
            S(:,:,i) = (W0(:,(i*params.numberOfSources/2) - (params.numberOfSources/2 -1):(i*params.numberOfSources/2)) * ...
                HTest{1,JNRIndex}((i*params.numberOfSources/2) - (params.numberOfSources/2 -1):(i*params.numberOfSources/2),:)./ (W0*HTest{1,JNRIndex})).*PxxTest{1,JNRIndex};
            
            xHat(:,i,JNRIndex,loopIndex) = istft(S(:,:,i), params.fs, 'Window', params.window, 'OverlapLength', params.overlap, 'FFTLength', params.nfft);
            
        end
    end
end

save(['.' filesep 'data' filesep 'nmf_testing_03.mat'], 'xHat', 'mixtureSignal', 'JNRVector', 'S');

rmpath(['.' filesep 'data']);
rmpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
rmpath(['..' filesep 'signalsGeneration' filesep]);
