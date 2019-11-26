clear;
clc;
close all;

addpath(['..' filesep 'Sigtools' filesep])
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
addpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['..' filesep 'Misc'])
addpath(['.' filesep 'data']);

load nmf_training_01.mat;
load sim_params_2.mat;

monteCarloLoops = 1;
SNR = -25;
alpha = 4.5e11;
deltaT = 12e-6;
params.JNRVector = [-10 0 10 20 30];

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
params.numberOfIterations = 1000;
params.tolChange = 1e-6;
params.tolError = 1e-6;
params.repetitions = 1;

t = -10e-6 - 1/params.fs:1/params.fs:20e-6 - 1/params.fs;
gaussFun = @(alpha, deltaT, t) exp((-alpha/2).*t.^2) + exp((-alpha/2).*(t - deltaT).^2);
dme =  gaussFun(alpha, deltaT, t).';

numberOfRawSamples = floor(paramsSignal.Freqsamp*paramsSignal.Intetime);
delay = 100e-6;
totalSamples = numberOfRawSamples;
numberOfZeros = numberOfRawSamples - length(dme);

for loopIndex = 1:monteCarloLoops
    loopIndex
    paramsSignal.numberOfGPSSignals = 1;
    paramsSignal.FreqDopp = 3e3;
    params.numberOfSources = params.numberOfSources*2;
    params.init = 'custom';
    params.W0 = W0;
    params.transform = false;
    
    %------------------Using different GPS and noise signals
    
    interferenceSignal = [zeros(numberOfZeros/4, 1);dme;zeros(numberOfZeros*3/4, 1)];
    Timeofthisloop = 0:totalSamples-1;
    Carrphase = mod(2*pi*(paramsSignal.FreqDopp)*Timeofthisloop/paramsSignal.Freqsamp,2*pi);
    Carrier = exp(1i*Carrphase).';
    
    interferenceSignal = interferenceSignal.*Carrier;
    paramsSignal.FreqDopp = 2e3;
    GPSSignals = GPSGen(paramsSignal);
    GPSSignals = GPSSignals(1:numberOfRawSamples,:);
    
    GPSSignals = [zeros(round(delay*params.fs), 1);GPSSignals];
    GPSSignalsPower = pow_eval(GPSSignals);
    interferenceSignal = [zeros(round(delay*params.fs), 1);interferenceSignal];
    interferenceSignalPower = pow_eval(interferenceSignal);

    noise = randn(length(GPSSignals), 1) + 1j*randn(length(GPSSignals), 1);
    noisePower = pow_eval(noise);
    mixtureSignal = zeros(length(noise), length(params.JNRVector));
    
    for JNRIndex = 1:length(params.JNRVector)
        GPSSignalsAux = GPSSignals;
        interferenceSignalAux = interferenceSignal;
        GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
        mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
        interferenceSignalAux = interferenceSignalAux*sqrt(noisePower*10^(params.JNRVector(JNRIndex)/10)/interferenceSignalPower);
        mixtureSignal(:,JNRIndex) = mixtureGPS + interferenceSignalAux;
    end
    
    [WTest, HTest, errorTest, PxxTest, f, t] = nmf_eval_v2(mixtureSignal, params);
    S = zeros([size(PxxTest{1,1}) 2]);
    xHat = zeros(length(mixtureSignal), 2, length(params.JNRVector));
    
    for JNRIndex = 1:length(params.JNRVector)
        JNRIndex
        for i = 1:2
            S(:,:,i) = (W0(:,(i*params.numberOfSources/2) - (params.numberOfSources/2 -1):(i*params.numberOfSources/2)) * ...
                HTest{1,JNRIndex}((i*params.numberOfSources/2) - (params.numberOfSources/2 -1):(i*params.numberOfSources/2),:)./ (W0*HTest{1,JNRIndex})).*PxxTest{1,JNRIndex};
            
            xHat(:,i,JNRIndex) = istft(S(:,:,i), params.fs, 'Window', params.window, 'OverlapLength', params.overlap, 'FFTLength', params.nfft);
            
        end
    end
end

save(['.' filesep 'data' filesep 'nmf_testing_03.mat'], 'xHat', 'mixtureSignal', 'WTest', 'HTest', 'errorTest', 'PxxTest', 'mixtureSignal');

rmpath(['.' filesep 'data']);
rmpath(['..' filesep 'Misc'])
rmpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);