clear;
clc;
close all;

addpath(['..' filesep 'Sigtools' filesep])
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
addpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])

load sim_params_L5_1.mat;
% load(['.' filesep 'DME signals' filesep 'signals_50MHz_snap_3.mat']);
matObjDME = matfile(['.' filesep 'DME signals' filesep 'signals_50MHz_snap_10_FL390.mat']);
matObjSOI = matfile(['.' filesep 'data' filesep 'W0_test_number_components_SOI_60_200.mat']);
matObjGPS = matfile(['.' filesep 'GPS signals' filesep 'GPS_L5.mat']);

SNR = -20;
params.JNRVector = 10:5:40;

numberOfComponentsVector = 60;
params.numberOfSources = numberOfComponentsVector;
numberOfSources = 2;

JNRVector = params.JNRVector;
params.fs = paramsSignal.Freqsamp;
params.nfft = 256;
params.nperseg = 256;
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
params.verbose = false;
params.alg = 'vanilla_semi';
params.numberOfComponentsPerSource = [params.numberOfSources params.numberOfSources];
params.numberOfSources = params.numberOfSources*numberOfSources;
params.W0 = matObjSOI.WSOI;
params.W0 = params.W0{1,1};

numPeriods = 1;
numberOfRawSamples = floor(paramsSignal.Freqsamp*paramsSignal.Intetime)*numPeriods;
delay = 10e-6;
totalSamples = numberOfRawSamples;

GPSSignalsPower = pow_eval(matObjGPS.GPSSignals);
interferenceSignalBuffer = buffer(matObjDME.DME, 1e-3*params.fs, 0, 'nodelay');
monteCarloLoops = 200;

mixtureSignal = zeros(totalSamples, length(params.JNRVector), monteCarloLoops);
xHat = zeros(totalSamples, length(params.JNRVector), length(numberOfComponentsVector), monteCarloLoops);

for loopIndex = 1:monteCarloLoops
    loopIndex
        
    noise = randn(paramsSignal.Intenumb, 1) + 1j*randn(paramsSignal.Intenumb, 1);
    noisePower = pow_eval(noise);    
    
    interferenceSignal = interferenceSignalBuffer(:,loopIndex);
    interferenceSignalPower = pow_eval(interferenceSignal);

    for JNRIndex = 1:length(params.JNRVector)
        GPSSignalsAux = matObjGPS.GPSSignals;
        interferenceSignalAux = interferenceSignal;
        GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
        mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
        interferenceSignalAux = interferenceSignalAux*sqrt(noisePower*10^(params.JNRVector(JNRIndex)/10)/interferenceSignalPower);
        mixtureSignal(:,JNRIndex,loopIndex) = mixtureGPS + interferenceSignalAux;
    end
    for numberOfComponentsIndex = 1:length(numberOfComponentsVector)
        
        [WTest, HTest, errorTest, PxxTest, f, t] = nmf_eval_v2(mixtureSignal(:,:,loopIndex), params);
        
        for JNRIndex = 1:length(params.JNRVector)
            S = (WTest{1,JNRIndex}(:,numberOfComponentsVector(numberOfComponentsIndex)+1:end) * ...
                HTest{1,JNRIndex}(numberOfComponentsVector(numberOfComponentsIndex)+1:end,:)./ (WTest{1,JNRIndex}*HTest{1,JNRIndex})).*PxxTest{1,JNRIndex};
            xHat(:,JNRIndex,numberOfComponentsIndex,loopIndex) = istft(S, params.fs, 'Window', params.window, 'OverlapLength', params.overlap, 'FFTLength', params.nfft);
        end
    end
end

xHat = single(xHat);
mixtureSignal = single(mixtureSignal);

save(['.' filesep 'data' filesep 'nmf_testing_12.mat'], 'xHat', 'mixtureSignal', 'JNRVector', 'numberOfComponentsVector');

rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
rmpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])