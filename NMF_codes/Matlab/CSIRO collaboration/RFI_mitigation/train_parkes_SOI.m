%This script trains the NMF using the new dataset for the 1152 channel
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
    
addpath(['..' filesep '..' filesep 'Sigtools' filesep 'NMF_algorithms']);

frequency = '1152';
fileName = '2018-09-01-09_52_55_0000000000000000.000000'; %test
year = '2018';
fileNumber = 1;

filePath2 = ['E:' filesep 'Parkes data' filesep year filesep frequency filesep 'signal frames' filesep];
matObj = matfile(['..' filesep 'Labels' filesep year filesep frequency filesep fileName '_labels.mat']);

falseLabels = find(~matObj.interferenceDetFlag);

numberOfSignalFrames = 20;
nonInterferenceFramesA = zeros(1280000, numberOfSignalFrames);
nonInterferenceFramesB = zeros(1280000, numberOfSignalFrames);

monteCarloLoops = 1;
params.JNRVector = 0;
params.fs = 128e6;
params.nfft = 256;
params.nperseg = 256;
params.overlap = params.nperseg/2;
params.window = ones(params.nperseg, 1);
params.hop_size = params.nperseg - params.overlap;
params.specType = 'power';
params.numberOfSources = 10;
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 500;
params.tolChange = 1e-6;
params.tolError = 1e-6;
params.repetitions = 1;
params.verbose = 1;
params.alg = 'vanilla';
params.centered = false;

for i = 1:numberOfSignalFrames
    load([filePath2 fileName '_' num2str(falseLabels(i)) '.mat']);
    nonInterferenceFramesA(:,i) = parkesSignalA;
    nonInterferenceFramesB(:,i) = parkesSignalB;
end

nonInterferenceFramesA = nonInterferenceFramesA(:);
nonInterferenceFramesB = nonInterferenceFramesB(:);

for loopIndex = 1:monteCarloLoops
    [WSignalA, ~, ~, ~, ~, ~] = nmf_eval_v2(nonInterferenceFramesA, params);  
    [WSignalB, ~, ~, ~, ~, ~] = nmf_eval_v2(nonInterferenceFramesB, params); 
    W0SOIA = [WSignalA{1,1}];
    W0SOIB = [WSignalB{1,1}];
end

save(['.' filesep 'data' filesep year filesep frequency filesep 'W_SOI_' num2str(fileNumber) '_'  fileName '.mat'],...
    'W0SOIA', 'W0SOIB', 'params', 'fileName', '-v7.3');             %saving trained matrix W

rmpath(['..' filesep '..' filesep 'Sigtools' filesep 'NMF_algorithms']);