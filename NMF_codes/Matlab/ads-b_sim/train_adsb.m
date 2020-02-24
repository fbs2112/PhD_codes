clear;
clc;
close all;
    
addpath(['..' filesep 'Sigtools' filesep])
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
addpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms']);
addpath(['.' filesep 'data']);

load adsb_signal;
load sim_params_2.mat;

monteCarloLoops = 1;

params.JNRVector = 0;

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
params.tolChange = 1e-3;
params.tolError = 1e-3;
params.repetitions = 1;
params.verbose = 1;

rcv = rcv(1:end/2);

totalSamples = length(rcv);
JNR = 20;
SNR = -20;

% numberOfZeros = numberOfRawSamples - length(dme);
% 
% interferenceSignal = [zeros(numberOfZeros/2, 1);dme;zeros(numberOfZeros/2, 1)];
% Timeofthisloop = 0:totalSamples-1;               
% Carrphase = mod(2*pi*(paramsSignal.FreqDopp)*Timeofthisloop/paramsSignal.Freqsamp,2*pi);
% Carrier = exp(1i*Carrphase).';
% 
% interferenceSignal = interferenceSignal.*Carrier;
% interferenceSignalPower = pow_eval(interferenceSignal);
% 
% paramsSignal.numberOfGPSSignals = 1; %overriding number of GPS signals from configuration file
% GPSSignals = GPSGen(paramsSignal);
% GPSSignals = GPSSignals(1:numberOfRawSamples,:);
% GPSSignalsPower = pow_eval(GPSSignals);

numberOfRawSamples = floor(paramsSignal.Freqsamp*paramsSignal.Intetime);
totalSamples = numberOfRawSamples;

paramsSignal.numberOfGPSSignals = 1; %overriding number of GPS signals from configuration file
GPSSignals = GPSGen(paramsSignal);
GPSSignals = GPSSignals(1:numberOfRawSamples,:);
GPSSignalsPower = pow_eval(GPSSignals);

for loopIndex = 1:monteCarloLoops
    loopIndex
    noise = single(randn(totalSamples, 1) + 1j*randn(totalSamples, 1));
    noisePower = pow_eval(noise);
    
    GPSSignalsAux = GPSSignals;
    GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
    mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;

    interferenceSignal = rcv;
    interferencePower = pow_eval(interferenceSignal);
    
    interferenceMultiplier = sqrt(noisePower*10.^(JNR/10)./interferencePower);
    mixtureRFI = interferenceSignal.*interferenceMultiplier;

    [WInterf, HInterf, errorInterfTrain, PxxInterf, ~, ~] = nmf_eval_v2(mixtureRFI, params);
    [WSignal, HSignal, errorSignalTrain, PxxSignal, ~, ~] = nmf_eval_v2(noise, params);
    
    W0 = [WInterf{1,1} WSignal{1,1}];
end

save(['.' filesep 'data' filesep 'nmf_training_01.mat'], 'W0');

rmpath(['.' filesep 'data']);
rmpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms']);
rmpath(['..' filesep 'Sigtools' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
rmpath(['..' filesep 'signalsGeneration' filesep]);