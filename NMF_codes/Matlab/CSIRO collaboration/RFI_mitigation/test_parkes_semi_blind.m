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

addpath(['..' filesep '..' filesep 'Sigtools' filesep])
addpath(['..' filesep '..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep '..' filesep 'Sigtools' filesep 'NMF_algorithms'])

frequency = '1152';
fileName = '2018-09-01-09_52_55_0000000000000000.000000'; %test
year = '2018';

load(['..' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep fileName '_labels.mat']);
load(['.' filesep 'data' filesep year filesep frequency filesep 'results1.mat']); %trained W0

params.numberOfComponentsPerSource = [10 10];
params.verbose = true;
params.alg = 'vanilla_semi';
params.mu = [0 0];
params.numberOfSources = sum(params.numberOfSources);

paramsA = params;
paramsB = params;
paramsA.W0 = W0A;
paramsB.W0 = W0B;

trueLabels = find(interferenceDetFlag);
trueLabels = sort(trueLabels);
trueLabels = trueLabels(trueLabels <= 100);
falseLabels = find(~interferenceDetFlag);
falseLabels = sort(falseLabels);
falseLabels = falseLabels(falseLabels < 100);

numberOfSignalFrames = length(trueLabels);

xHatA = zeros(1280000, 2, length(params.JNRVector), 100);
xHatB = zeros(1280000, 2, length(params.JNRVector), 100);

WA = cell(1,numberOfSignalFrames);
HA = cell(1,numberOfSignalFrames);
WB = cell(1,numberOfSignalFrames);
HB = cell(1,numberOfSignalFrames);

for loopIndex = 1:numberOfSignalFrames
    disp(['Loop index: ' num2str(loopIndex)]);
    load(['..' filesep 'data' filesep 'dataParkesNew' filesep frequency filesep fileName '_' num2str(trueLabels(loopIndex)) '.mat']);
    interferenceFramesA = parkesSignalA;
    interferenceFramesB = parkesSignalB;
    [WTestA, HTestA, errorCell, PxxA, f, t] = nmf_eval_v2(interferenceFramesA, paramsA);
    [WTestB, HTestB, ~, PxxB, ~, ~] = nmf_eval_v2(interferenceFramesB, paramsB);
%     aux = WTest{1,1}(:,1:10);
%     aux(aux < max(max(aux))/2) = 0;
    SA = zeros([size(PxxA{1,1}) 2]);
    WA{1,loopIndex} = WTestA;
    HA{1,loopIndex} = HTestA;
%     WTest{1,1}(:,1:10) = aux;

    %Polarisaton A---------------------------------------------------------
    for i = 1:2
        %------Wiener solution to reconstruct the complex spectrograms
        SA(:,:,i) = (WTestA{1,1}(:,((i-1)*params.numberOfSources/2) + 1:(i*params.numberOfSources/2)) * ...
            HTestA{1,1}(((i-1)*params.numberOfSources/2) + 1:(i*params.numberOfSources/2),:)./ (WTestA{1,1}*HTestA{1,1})).*PxxA{1,1};
        %--------------------------------------------------------------
        %Time domain reconstruction via iSTFT--------------------------
        xHatA(:,i,1,trueLabels(loopIndex)) = istft(SA(:,:,i), params.fs, 'Window', ....
            params.window, 'OverlapLength', params.overlap, 'FFTLength', params.nfft, 'centered', false);
        %--------------------------------------------------------------
    end
    
    %----------------------------------------------------------------------
    
    %Polarisaton B---------------------------------------------------------
    SB = zeros([size(PxxB{1,1}) 2]);
    WB{1,loopIndex} = WTestB;
    HB{1,loopIndex} = HTestB;
    for i = 1:2
        %------Wiener solution to reconstruct the complex spectrograms
        SB(:,:,i) = (WTestB{1,1}(:,((i-1)*params.numberOfSources/2) + 1:(i*params.numberOfSources/2)) * ...
            HTestB{1,1}(((i-1)*params.numberOfSources/2) + 1:(i*params.numberOfSources/2),:)./ (WTestB{1,1}*HTestB{1,1})).*PxxB{1,1};
        %--------------------------------------------------------------
        %Time domain reconstruction via iSTFT--------------------------
        xHatB(:,i,1,trueLabels(loopIndex)) = istft(SB(:,:,i), params.fs, 'Window',...
            params.window, 'OverlapLength', params.overlap, 'FFTLength', params.nfft, 'centered', false);
        %--------------------------------------------------------------
    end
    %----------------------------------------------------------------------
end

save(['.' filesep 'data' filesep year filesep frequency filesep 'results2.mat'], 'xHatA', 'xHatB', 'WA', 'HA', 'WB', 'HB');

rmpath(['..' filesep '..' filesep 'Sigtools' filesep])
rmpath(['..' filesep '..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep '..' filesep 'Sigtools' filesep 'NMF_algorithms'])