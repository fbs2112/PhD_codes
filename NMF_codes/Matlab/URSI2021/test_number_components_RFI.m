clear;
clc;
close all;
    
addpath(['..' filesep 'Sigtools' filesep])
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
addpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['..' filesep 'Third party' filesep 'Dani Pascual' filesep 'source']);
addpath(['..' filesep 'Third party' filesep 'Dani Pascual' filesep 'source' filesep 'functions' filesep 'signal_processing']);
addpath(['..' filesep 'Third party' filesep 'Dani Pascual' filesep 'prn_codes']);

load sim_params_L5_1.mat;
load(['.' filesep 'DME signals' filesep 'signals_50MHz_snap_3_FL390.mat']);
load(['.' filesep 'GPS signals' filesep 'GPS_L5_train.mat']);

monteCarloLoops = 1;
SNR = -20;

params.JNRVector = 0;
params.fs = paramsSignal.Freqsamp;
params.nfft = 1024;
params.nperseg = 1024;
params.overlap = params.nperseg - 1;
params.hop_size = params.nperseg - params.overlap;
params.window = ones(params.nperseg, 1);
params.specType = 'power';
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 500;
params.tolChange = 1e-6;
params.tolError = 1e-6;
params.repetitions = 1;
params.verbose = 1;

numPeriods = 2;

numberOfRawSamples = floor(paramsSignal.Freqsamp*paramsSignal.Intetime)*numPeriods;
totalSamples = numberOfRawSamples;

interferenceSignal = DME(1:totalSamples).';
interferenceSignalPower = pow_eval(interferenceSignal);
numberOfComponentsVec = 400;

for loopIndex = 1:monteCarloLoops
    for i = 1:length(numberOfComponentsVec)
        params.numberOfSources = numberOfComponentsVec(i);
        [WInterf, HInterf, errorInterfTrain, PxxInterf, f, t] = nmf_eval_v2(interferenceSignal, params);
        WRFI{i,1} = WInterf{1,1};
    end
end
save(['.' filesep 'data' filesep 'W0_test_number_components_RFI_400_1024_FL390.mat'], 'WRFI');

rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
rmpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep 'Third party' filesep 'Dani Pascual' filesep 'source']);
rmpath(['..' filesep 'Third party' filesep 'Dani Pascual' filesep 'source' filesep 'functions' filesep 'signal_processing']);
rmpath(['..' filesep 'Third party' filesep 'Dani Pascual' filesep 'prn_codes']);