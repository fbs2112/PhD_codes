%This script trains the NMF according to a given ADS-B signal and the
%signal of interest (CW)
%
%This code was created by Felipe Barboza da Silva
% Copyright (c) School of Engineering. Macquarie University - Australia.
% All rights reserved. 2020.
% DISCLAIMER:
%     This code is strictly private, confidential and personal to its recipients
%     and should not be copied, distributed or reproduced in whole or in part,
%     nor passed to any third party.
% ** Do not duplicate or distribute without written permission from the owners. **

clear;
clc;
close all;
    
addpath(['..' filesep 'Sigtools' filesep])
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms']);
addpath(['.' filesep 'data']);

load adsb_signal;

monteCarloLoops = 1;
params.JNRVector = 0;
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
params.verbose = 1;
JNR = 20;

rcv = rcv(1:end/2);
rcv = rcv(1.5e4:1.75e4);
totalSamples = length(rcv);

bandwidthVector = 0;
initialFrequency = 1e6;

for loopIndex = 1:monteCarloLoops
    paramsSignal.Freqsamp = params.fs;    
    paramsSignal.Intenumb = length(rcv);
    paramsSignal.Noneperiod = paramsSignal.Intenumb;                           % number of samples with a sweep time
    paramsSignal.IFmin = initialFrequency;                                     % start frequency
    paramsSignal.IFmax = bandwidthVector + initialFrequency;                   % end frequency
    paramsSignal.foneperiod(1:paramsSignal.Noneperiod) = linspace(paramsSignal.IFmin, paramsSignal.IFmax, paramsSignal.Noneperiod);
    paramsSignal.Initphase = 0;
    
    signalOfInterest = interferenceGen(paramsSignal);                      %CW generation
    signalOfInterestPower = pow_eval(signalOfInterest);
    interferenceSignal = rcv;
    interferencePower = pow_eval(interferenceSignal);
    
    interferenceMultiplier = sqrt(signalOfInterestPower*10.^(JNR/10)./interferencePower);
    mixtureRFI = interferenceSignal.*interferenceMultiplier;

    [WInterf, HInterf, errorInterfTrain, PxxInterf, ~, ~] = nmf_eval_v2(mixtureRFI, params);                %ADS-B frequency content evaluation
    [WSignal, HSignal, errorSignalTrain, PxxSignal, ~, ~] = nmf_eval_v2(single(signalOfInterest), params);  %CW frequency content evaluation
    
    W0 = [WInterf{1,1} WSignal{1,1}];
end

save(['.' filesep 'data' filesep 'nmf_training_03.mat'], 'W0');             %saving trained matrix W

rmpath(['.' filesep 'data']);
rmpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms']);
rmpath(['..' filesep 'Sigtools' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep]);