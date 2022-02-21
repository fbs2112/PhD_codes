clear;
clc;
close all;

load(['.' filesep 'data' filesep 'nmf_testing_13_1']);
xOut(:,1,1,1,:) = xHat;
mixtureSignalMerged(:,1,:) = mixtureSignal;

load(['.' filesep 'data' filesep 'nmf_testing_13_2']);
xOut(:,1,2,1,:) = xHat;
mixtureSignalMerged(:,2,:) = mixtureSignal;

load(['.' filesep 'data' filesep 'nmf_testing_13_3']);
xOut(:,1,3,1,:) = xHat;
mixtureSignalMerged(:,3,:) = mixtureSignal;

load(['.' filesep 'data' filesep 'nmf_testing_13_4']);
xOut(:,1,4,1,:) = xHat;
mixtureSignalMerged(:,4,:) = mixtureSignal;

load(['.' filesep 'data' filesep 'nmf_testing_13_5']);
xOut(:,1,5,1,:) = xHat;
mixtureSignalMerged(:,5,:) = mixtureSignal;

load(['.' filesep 'data' filesep 'nmf_testing_13_6']);
xOut(:,1,6,1,:) = xHat;
mixtureSignalMerged(:,6,:) = mixtureSignal;

load(['.' filesep 'data' filesep 'nmf_testing_13_7']);
xOut(:,1,7,1,:) = xHat;
mixtureSignalMerged(:,7,:) = mixtureSignal;

JNRVector = [10 15 20 25 30 35 40];

save(['.' filesep 'data' filesep 'nmf_testing_13.mat'], 'xOut', 'mixtureSignalMerged', 'JNRVector', 'numberOfComponentsVector');
