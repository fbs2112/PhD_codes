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

addpath(['..' filesep '..' filesep 'Sigtools' filesep 'NMF_algorithms'])

frequency = '1152';
fileName = '2018-09-01-09_52_55_0000000000000000.000000'; %test
year = '2018';
filePath2 = ['F:' filesep 'Parkes data' filesep year filesep frequency filesep 'signal frames' filesep];

fileNumberSOI = 1;
fileNumberRFI = 2;

matObj = matfile(['..' filesep 'Labels' filesep year filesep frequency filesep fileName '_labels.mat']);
trueLabels = find(matObj.interferenceDetFlag == 1 | matObj.interferenceDetFlag == 3 | matObj.interferenceDetFlag == 5 | matObj.interferenceDetFlag == 6);

matObjSOI = matfile(['.' filesep 'data' filesep year filesep frequency filesep 'W_SOI_' num2str(fileNumberSOI) '_'  fileName '.mat']);

W0SOIA = matObjSOI.W0SOIA;
W0SOIB = matObjSOI.W0SOIB;

params = matObjSOI.params;

W0RFIA = zeros(params.nfft, 1);
W0RFIB = W0RFIA;

W0RFIA(133-7:133+7,:) = 1;
W0RFIB(133-7:133+7,:) = 1;

params.numberOfComponentsPerSource = [1 10];
params.verbose = false;
params.numberOfSources = sum(params.numberOfComponentsPerSource);

params.init = 'custom';
params.transform = false;
paramsA = params;
paramsB = params;
paramsA.W0 = [W0RFIA W0SOIA];
paramsB.W0 = [W0RFIB W0SOIB];

numberOfSignalFrames = find(trueLabels >= 100, 1, 'first');

xHatA = zeros(1280000, 2, length(params.JNRVector), 100);
xHatB = zeros(1280000, 2, length(params.JNRVector), 100);

WA = cell(1,numberOfSignalFrames);
HA = cell(1,numberOfSignalFrames);
WB = cell(1,numberOfSignalFrames);
HB = cell(1,numberOfSignalFrames);

for loopIndex = 1:numberOfSignalFrames
    disp(['Loop index: ' num2str(loopIndex)]);
    load([filePath2 fileName '_' num2str(trueLabels(loopIndex)) '.mat']);
    interferenceFramesA = parkesSignalA;
    interferenceFramesB = parkesSignalB;
    
    interferenceFramesABuffer = buffer(interferenceFramesA, length(interferenceFramesA)/1, 0, 'nodelay');
    interferenceFramesBBuffer = buffer(interferenceFramesB, length(interferenceFramesB)/1, 0, 'nodelay');
    
    xHatAAux = zeros(size(interferenceFramesABuffer, 1), size(interferenceFramesABuffer, 2), 2);
    xHatBAux = zeros(size(interferenceFramesABuffer, 1), size(interferenceFramesABuffer, 2), 2);
    for signalIndex = 1:size(interferenceFramesABuffer, 2)
        [WTestA, HTestA, errorCell, PxxA, f, t] = nmf_eval_v2(interferenceFramesABuffer(:,signalIndex), paramsA);
        [WTestB, HTestB, ~, PxxB, ~, ~] = nmf_eval_v2(interferenceFramesBBuffer(:,signalIndex), paramsB);
        
        %Polarisation A---------------------------------------------------------
        SA = zeros([size(PxxA{1,1}) 2]);
        WA{1,loopIndex} = WTestA;
        HA{1,loopIndex} = HTestA;
        %------Wiener solution to reconstruct the complex spectrograms
        SA(:,:,1) = (WTestA{1,1}(:,1:params.numberOfComponentsPerSource(1)) * ...
            HTestA{1,1}(1:params.numberOfComponentsPerSource(1),:) ./ (WTestA{1,1}*HTestA{1,1})).*PxxA{1,1};
        SA(:,:,2) = (WTestA{1,1}(:,params.numberOfComponentsPerSource(1)+1:end) * ...
            HTestA{1,1}(params.numberOfComponentsPerSource(1)+1:end,:) ./ (WTestA{1,1}*HTestA{1,1})).*PxxA{1,1};
        %--------------------------------------------------------------
        for i = 1:2
       
            %Time domain reconstruction via iSTFT--------------------------
             xHatAAux(:,signalIndex,i) = istft(SA(:,:,i), params.fs, 'Window', ....
                params.window, 'OverlapLength', params.overlap, 'FFTLength', params.nfft, 'centered', false);
            %--------------------------------------------------------------
        end
        
        %----------------------------------------------------------------------
        
        %Polarisation B---------------------------------------------------------
        SB = zeros([size(PxxB{1,1}) 2]);
        WB{1,loopIndex} = WTestB;
        HB{1,loopIndex} = HTestB;
        %------Wiener solution to reconstruct the complex spectrograms
        SB(:,:,1) = (WTestB{1,1}(:,1:params.numberOfComponentsPerSource(1)) * ...
            HTestB{1,1}(1:params.numberOfComponentsPerSource(1),:) ./ (WTestB{1,1}*HTestB{1,1})).*PxxB{1,1};
        SB(:,:,2) = (WTestB{1,1}(:,params.numberOfComponentsPerSource(1)+1:end) * ...
            HTestB{1,1}(params.numberOfComponentsPerSource(1)+1:end,:) ./ (WTestB{1,1}*HTestB{1,1})).*PxxB{1,1};
        %--------------------------------------------------------------
        for i = 1:2
            
            %Time domain reconstruction via iSTFT--------------------------
             xHatBAux(:,signalIndex,i)  = istft(SB(:,:,i), params.fs, 'Window',...
                params.window, 'OverlapLength', params.overlap, 'FFTLength', params.nfft, 'centered', false);
            %--------------------------------------------------------------
        end
    end
    
    xHatA(:,:,1,trueLabels(loopIndex)) = reshape(xHatAAux, [], 2);
    xHatB(:,:,1,trueLabels(loopIndex)) = reshape(xHatBAux, [], 2);
end

save(['.' filesep 'data' filesep year filesep frequency filesep 'results9.mat'], 'xHatA', 'xHatB', '-v7.3');

rmpath(['..' filesep '..' filesep 'Sigtools' filesep 'NMF_algorithms'])