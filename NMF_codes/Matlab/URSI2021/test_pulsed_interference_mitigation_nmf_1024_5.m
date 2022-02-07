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
% load(['.' filesep 'DME signals' filesep 'signals_50MHz_snap_3.mat']);
load(['.' filesep 'DME signals' filesep 'signals_50MHz_snap_10_FL390.mat']);
load(['.' filesep 'data' filesep 'W0_test_number_components_RFI_400_1024_FL390.mat']);
load(['.' filesep 'data' filesep 'W0_test_number_components_SOI_400_1024.mat']);
load(['.' filesep 'GPS signals' filesep 'GPS_L5.mat']);

SNR = -20;
params.JNRVector = [30];

JNRVector = params.JNRVector;
params.fs = paramsSignal.Freqsamp;
params.nfft = 1024;
params.nperseg = 1024;
params.overlap = params.nperseg - 1;
params.hop_size = params.nperseg - params.overlap;
params.window = ones(params.nperseg, 1);
params.specType = 'power';
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 300;
params.tolChange = 1e-6;
params.tolError = 1e-6;
params.repetitions = 1;
params.init = 'custom';
params.transform = false;
params.verbose = false;

numberOfComponentsVector = 400;
numPeriods = 1;
numberOfRawSamples = floor(paramsSignal.Freqsamp*paramsSignal.Intetime)*numPeriods;
delay = 10e-6;
totalSamples = numberOfRawSamples;

GPSSignalsPower = pow_eval(GPSSignals);
interferenceSignalBuffer = buffer(DME, 1e-3*params.fs, 0, 'nodelay');
monteCarloLoops = 200;

mixtureSignal = zeros(totalSamples, length(params.JNRVector), monteCarloLoops);
xHat = zeros(totalSamples, 2, length(params.JNRVector), length(numberOfComponentsVector), monteCarloLoops);

for loopIndex = 1:monteCarloLoops
    loopIndex
        
    noise = randn(length(GPSSignals), 1) + 1j*randn(length(GPSSignals), 1);
    noisePower = pow_eval(noise);    
    
    interferenceSignal = interferenceSignalBuffer(:,loopIndex);
    interferenceSignalPower = pow_eval(interferenceSignal);

    for JNRIndex = 1:length(params.JNRVector)
        GPSSignalsAux = GPSSignals;
        interferenceSignalAux = interferenceSignal;
        GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
        mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
        interferenceSignalAux = interferenceSignalAux*sqrt(noisePower*10^(params.JNRVector(JNRIndex)/10)/interferenceSignalPower);
        mixtureSignal(:,JNRIndex,loopIndex) = mixtureGPS + interferenceSignalAux;
    end
    for numberOfComponentsIndex = 1:length(numberOfComponentsVector)
        
        params.W0 = [WRFI{numberOfComponentsIndex,1} WSOI{numberOfComponentsIndex,1}];
        params.numberOfSources = size(params.W0, 2);
        [WTest, HTest, errorTest, PxxTest, f, t] = nmf_eval_v2(mixtureSignal(:,:,loopIndex), params);
        S = zeros([size(PxxTest{1,1}) 2]);
        
        for JNRIndex = 1:length(params.JNRVector)
            S(:,:,1) = (params.W0(:,1:numberOfComponentsVector(numberOfComponentsIndex)) * ...
                HTest{1,JNRIndex}(1:numberOfComponentsVector(numberOfComponentsIndex),:)./ (params.W0*HTest{1,JNRIndex})).*PxxTest{1,JNRIndex};
            
            S(:,:,2) = (params.W0(:,numberOfComponentsVector(numberOfComponentsIndex)+1:end) * ...
                HTest{1,JNRIndex}(numberOfComponentsVector(numberOfComponentsIndex)+1:end,:)./ (params.W0*HTest{1,JNRIndex})).*PxxTest{1,JNRIndex};
            for i = 1:2
                xHat(:,i,JNRIndex,numberOfComponentsIndex,loopIndex) = istft(S(:,:,i), params.fs, 'Window', params.window, 'OverlapLength', params.overlap, 'FFTLength', params.nfft);
            end
        end
    end
end

xHat = single(xHat(:,2,:,:,:));
mixtureSignal = single(mixtureSignal);

save(['.' filesep 'data' filesep 'nmf_testing_11_5.mat'], 'xHat', 'mixtureSignal', 'JNRVector', 'numberOfComponentsVector');

rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
rmpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep 'Third party' filesep 'Dani Pascual' filesep 'source']);
rmpath(['..' filesep 'Third party' filesep 'Dani Pascual' filesep 'source' filesep 'functions' filesep 'signal_processing']);
rmpath(['..' filesep 'Third party' filesep 'Dani Pascual' filesep 'prn_codes']);
