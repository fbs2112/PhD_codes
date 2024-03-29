clear;
clc;
close all;
    
addpath(['..' filesep 'Sigtools' filesep])
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
addpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])

load sim_params_2.mat;

monteCarloLoops = 1;
SNR = -25;
alpha = 4.5e11;
deltaT = 12e-6;
nbits = [2 4];

params.JNRVector = 0;
params.fs = paramsSignal.Freqsamp;
params.nfft = 256;
params.nperseg = 256;
params.overlap = params.nperseg - 1;
params.hop_size = params.nperseg - params.overlap;
params.window = ones(params.nperseg, 1);
params.specType = 'power';
params.numberOfSources = 5;
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 500;
params.tolChange = 1e-3;
params.tolError = 1e-3;
params.repetitions = 1;
params.verbose = 1;

t = -10e-6 - 1/params.fs:1/params.fs:20e-6 - 1/params.fs;
gaussFun = @(alpha, deltaT, t) exp((-alpha/2).*t.^2) + exp((-alpha/2).*(t - deltaT).^2);
dme = 3*gaussFun(alpha, deltaT, t).';

numberOfRawSamples = floor(paramsSignal.Freqsamp*paramsSignal.Intetime);
totalSamples = numberOfRawSamples;
numberOfZeros = numberOfRawSamples - length(dme);

interferenceSignal = [zeros(numberOfZeros/2, 1);dme;zeros(numberOfZeros/2, 1)];
% Timeofthisloop = 0:totalSamples-1;               
% Carrphase = mod(2*pi*(paramsSignal.FreqDopp)*Timeofthisloop/paramsSignal.Freqsamp,2*pi);
% Carrier = exp(1i*Carrphase).';

% interferenceSignal = interferenceSignal.*Carrier;
interferenceSignalPower = pow_eval(interferenceSignal);

paramsSignal.numberOfGPSSignals = 1; %overriding number of GPS signals from configuration file
paramsSignal.params.FreqDopp = 0;
GPSSignals = GPSGen(paramsSignal);
GPSSignals = GPSSignals(1:numberOfRawSamples,:);
GPSSignalsPower = pow_eval(GPSSignals);

W0 = zeros(params.nfft, params.numberOfSources*2, length(nbits));

for loopIndex = 1:monteCarloLoops
%     noise = randn(totalSamples, 1) + 1j*randn(totalSamples, 1);
    noise = randn(totalSamples, 1);

    noisePower = pow_eval(noise);

    GPSSignalsAux = GPSSignals;
    GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
    mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
    
    for nbitsIndex = 1:length(nbits)
        nbitsIndex
        interferenceSignalQuant = quantise_gps(interferenceSignal, nbits(nbitsIndex), 1);
        mixtureGPSQuant = quantise_gps(mixtureGPS, nbits(nbitsIndex), 1);

        [WInterf, HInterf, errorInterfTrain, PxxInterf, ~, ~] = nmf_eval_v2(interferenceSignalQuant, params);
        [WSignal, HSignal, errorSignalTrain, PxxSignal, ~, ~] = nmf_eval_v2(mixtureGPSQuant, params);
   
        W0(:,:,nbitsIndex) = [WInterf{1,1} WSignal{1,1}];
    end
end

save(['.' filesep 'data' filesep 'nmf_training_10.mat'], 'W0');

rmpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);