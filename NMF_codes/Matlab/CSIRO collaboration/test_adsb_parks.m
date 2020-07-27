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

load nmf_training_ADSB_01.mat;                                                  %reading training file
load ADSB_label.mat;

params.JNRVector = 10;                                                     %defining the jammer-to-noise ratio

params.fs = 128e6;
params.nfft = 256;
params.nperseg = 256;
params.overlap = params.nperseg/2;
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

%-----------------------Supervised NMF-------------------------------------
params.numberOfSources = params.numberOfSources*2;
params.init = 'custom';
params.W0 = W0;
params.transform = false;
params.verbose = true;
%--------------------------------------------------------------------------

trueLabels = find(interferenceDetFlag);
falseLabels = find(~interferenceDetFlag);

numberOfSignalFrames = 10;
monteCarloLoops = numberOfSignalFrames;

interferenceFrames = zeros(1280000, numberOfSignalFrames);
nonInterferenceFrames = zeros(1280000, numberOfSignalFrames);

for i = 1:numberOfSignalFrames
    load(['.' filesep 'data' filesep 'dataParks_' num2str(trueLabels(i + numberOfSignalFrames)) '.mat']);
    interferenceFrames(:,i) = parksSignal;
end

xHat = zeros(size(interferenceFrames, 1), 2, length(params.JNRVector), monteCarloLoops);

for loopIndex = 1:monteCarloLoops
        
    [WTest, HTest, errorTest, PxxTest, f, t] = nmf_eval_v2(interferenceFrames(:,loopIndex), params);  %NMF mitigation
    S = zeros([size(PxxTest{1,1}) 2]);
    
        for i = 1:2
            %------Wiener solution to reconstruct the complex spectrograms
            S(:,:,i) = (W0(:,(i*params.numberOfSources/2) - (params.numberOfSources/2 -1):(i*params.numberOfSources/2)) * ...
                HTest{1,1}((i*params.numberOfSources/2) - (params.numberOfSources/2 -1):(i*params.numberOfSources/2),:)./ (W0*HTest{1,1})).*PxxTest{1,1};
            %--------------------------------------------------------------
            %Time domain reconstruction via iSTFT--------------------------
            xHat(:,i,1,loopIndex) = istft(S(:,:,i), params.fs, 'Window', params.window, 'OverlapLength', params.overlap, 'FFTLength', params.nfft);
            %--------------------------------------------------------------
        end
    
end

save(['.' filesep 'data' filesep 'nmf_testing_ADSB_01.mat'], 'xHat');

rmpath(['.' filesep 'data']);
rmpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep]);