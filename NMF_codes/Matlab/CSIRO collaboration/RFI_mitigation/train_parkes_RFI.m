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
fileNumber = 2;

filePath2 = ['E:' filesep 'Parkes data' filesep year filesep frequency filesep 'signal frames' filesep];
matObj = matfile(['..' filesep 'Labels' filesep year filesep frequency filesep fileName '_labels.mat']);

trueLabels = find(matObj.interferenceDetFlag == 1 | matObj.interferenceDetFlag == 3 | matObj.interferenceDetFlag == 5 | matObj.interferenceDetFlag == 6);

M = readmatrix(['..' filesep 'Labels' filesep year filesep frequency filesep fileName '_labels_adsb_occurrences.csv']);

numberOfSignalFrames = 20;

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

interferenceFramesAAux2 = cell(numberOfSignalFrames, 1);
interferenceFramesBAux2 = cell(numberOfSignalFrames, 1);

for i = 1:numberOfSignalFrames
    load([filePath2 fileName '_' num2str(trueLabels(i)) '.mat']);
    aux = M(i,:);
    aux = aux(~isnan(aux));
    interferenceFramesAAux = cell(size(aux, 2)/2, 1);
    interferenceFramesBAux = cell(size(aux, 2)/2, 1);
    
    for j = 1:size(aux, 2)/2
        idxAux1 = round(aux((j-1)*2+1)*1e-3*params.fs);
        idxAux2 = round(aux(j*2)*1e-3*params.fs);
        interferenceFramesAAux{j} = parkesSignalA(idxAux1:idxAux2, 1);
        interferenceFramesBAux{j} = parkesSignalB(idxAux1:idxAux2, 1);
    end
    
    interferenceFramesAAux2{i} = vertcat(interferenceFramesAAux{:});
    interferenceFramesBAux2{i} = vertcat(interferenceFramesBAux{:});
end

interferenceFramesA = vertcat(interferenceFramesAAux2{:});
interferenceFramesB = vertcat(interferenceFramesBAux2{:});

for loopIndex = 1:monteCarloLoops
    [WSignalA, ~, ~, ~, ~, ~] = nmf_eval_v2(interferenceFramesA, params);  
    [WSignalB, ~, ~, ~, ~, ~] = nmf_eval_v2(interferenceFramesB, params); 
    W0RFIA = [WSignalA{1,1}];
    W0RFIB = [WSignalB{1,1}];
end

save(['.' filesep 'data' filesep year filesep frequency filesep 'W_RFI_' num2str(fileNumber) '_'  fileName '.mat'],...
    'W0RFIA', 'W0RFIB', 'params', 'fileName', '-v7.3');             %saving trained matrix W

rmpath(['..' filesep '..' filesep 'Sigtools' filesep 'NMF_algorithms']);