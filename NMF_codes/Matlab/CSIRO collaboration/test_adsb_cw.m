%This script uses the trained W matrix to mitigate interference, saving the
%reconstructed signal after the NMF algorithm, the corrupted received
%signal and the complex reconstructed spectrograms 
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
addpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['.' filesep 'data']);

load nmf_training_03.mat;                                                  %reading training file
load adsb_signal.mat;

monteCarloLoops = 1;
params.JNRVector = 10;                                                     %defining the jammer-to-noise ratio

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

rcv = rcv(end/2+1:end);
totalSamples = length(rcv);
interferencePower = pow_eval(rcv);

bandwidthVector = 0;
initialFrequency = 1e6;

%-----------------------Supervised NMF-------------------------------------
params.numberOfSources = params.numberOfSources*2;
params.init = 'custom';
params.W0 = W0;
params.transform = false;
params.verbose = true;
%--------------------------------------------------------------------------

mixtureSignal = single(zeros(totalSamples, length(params.JNRVector), monteCarloLoops));
xHat = zeros(length(mixtureSignal), 2, length(params.JNRVector), monteCarloLoops);

for loopIndex = 1:monteCarloLoops

    paramsSignal.Freqsamp = params.fs;
    paramsSignal.Intenumb = length(rcv);
    paramsSignal.Noneperiod = paramsSignal.Intenumb;                       % number of samples with a sweep time
    paramsSignal.IFmin = initialFrequency;                                 % start frequency
    paramsSignal.IFmax = bandwidthVector + initialFrequency;               % end frequency
    paramsSignal.foneperiod(1:paramsSignal.Noneperiod) = ...
        linspace(paramsSignal.IFmin, paramsSignal.IFmax, paramsSignal.Noneperiod);
    paramsSignal.Initphase = 0;
    
    signalOfInterest = interferenceGen(paramsSignal);
    signalOfInterest = signalOfInterest(1:length(rcv));
    signalOfInterestPower = pow_eval(signalOfInterest);
    
    for JNRIndex = 1:length(params.JNRVector)
        interferenceSignalAux = rcv;
        interferenceMultiplier = sqrt(signalOfInterestPower*10.^(params.JNRVector(JNRIndex)/10)./interferencePower);
        mixtureSignal(:,JNRIndex,loopIndex) = (interferenceSignalAux.*interferenceMultiplier) + single(signalOfInterest); 
    end
    
    [WTest, HTest, errorTest, PxxTest, f, t] = nmf_eval_v2(mixtureSignal(:,:,loopIndex), params);  %NMF mitigation
    S = zeros([size(PxxTest{1,1}) 2]);
    
    for JNRIndex = 1:length(params.JNRVector)
        for i = 1:2
            %------Wiener solution to reconstruct the complex spectrograms
            S(:,:,i) = (W0(:,(i*params.numberOfSources/2) - (params.numberOfSources/2 -1):(i*params.numberOfSources/2)) * ...
                HTest{1,JNRIndex}((i*params.numberOfSources/2) - (params.numberOfSources/2 -1):(i*params.numberOfSources/2),:)./ (W0*HTest{1,JNRIndex})).*PxxTest{1,JNRIndex};
            %--------------------------------------------------------------
            %Time domain reconstruction via iSTFT--------------------------
            xHat(:,i,JNRIndex,loopIndex) = istft(S(:,:,i), params.fs, 'Window', params.window, 'OverlapLength', params.overlap, 'FFTLength', params.nfft);
            %--------------------------------------------------------------
        end
    end
end

save(['.' filesep 'data' filesep 'nmf_testing_03.mat'], 'xHat', 'mixtureSignal', 'S');

rmpath(['.' filesep 'data']);
rmpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep]);