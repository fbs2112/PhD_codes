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
params.tolChange = 1e-6;
params.tolError = 1e-6;
params.repetitions = 1;
params.verbose = 1;

fs = 2.4e6;
rcv = rcv(1:end/2);
rcv = rcv(1.5e4:1.75e4);

totalSamples = length(rcv);
JNR = 20;
SNR = -20;

bandwidthVector = 0;
periodVector = 8.62e-6;
initialFrequency = 1e6;

for loopIndex = 1:monteCarloLoops
    loopIndex
    paramsSignal.Freqsamp = fs;    
    paramsSignal.Noneperiod = round(periodVector*fs);                   % number of samples with a sweep time
    paramsSignal.IFmin = initialFrequency;                                                  % start frequency
    paramsSignal.IFmax = bandwidthVector + initialFrequency;                    % end frequency
    paramsSignal.foneperiod(1:paramsSignal.Noneperiod) = linspace(paramsSignal.IFmin, paramsSignal.IFmax, paramsSignal.Noneperiod);
    paramsSignal.Initphase = 0;
    
    signalOfInterest = interferenceGen(paramsSignal);
    signalOfInterest = signalOfInterest(1:length(rcv));
    signalOfInterestPower = pow_eval(signalOfInterest);
%     GPSSignalsAux = GPSSignals;
%     GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
%     mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;

    interferenceSignal = rcv;
    interferencePower = pow_eval(interferenceSignal);
    
    interferenceMultiplier = sqrt(signalOfInterestPower*10.^(JNR/10)./interferencePower);
    mixtureRFI = interferenceSignal.*interferenceMultiplier;

    [WInterf, HInterf, errorInterfTrain, PxxInterf, ~, ~] = nmf_eval_v2(mixtureRFI, params);
    [WSignal, HSignal, errorSignalTrain, PxxSignal, ~, ~] = nmf_eval_v2(single(signalOfInterest), params);
    
    W0 = [WInterf{1,1} WSignal{1,1}];
end

save(['.' filesep 'data' filesep 'nmf_training_03.mat'], 'W0');

rmpath(['.' filesep 'data']);
rmpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms']);
rmpath(['..' filesep 'Sigtools' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
rmpath(['..' filesep 'signalsGeneration' filesep]);