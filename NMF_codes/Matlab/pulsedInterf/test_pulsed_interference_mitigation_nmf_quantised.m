clear;
clc;
close all;

addpath(['..' filesep 'Sigtools' filesep])
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
addpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['.' filesep 'data']);

load nmf_training_10.mat;
load sim_params_2.mat;

monteCarloLoops = 100;
SNR = -25;
nbits = [2 4];
alpha = 4.5e11;
deltaT = 12e-6;
params.JNRVector = [-5 0 10 30 50];
% params.JNRVector = [-5];

JNRVector = params.JNRVector; 
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
params.numberOfSources = params.numberOfSources*2;
params.init = 'custom';
params.transform = false;
params.verbose = false;

t = -10e-6 - 1/params.fs:1/params.fs:20e-6 - 1/params.fs;
gaussFun = @(alpha, deltaT, t) exp((-alpha/2).*t.^2) + exp((-alpha/2).*(t - deltaT).^2);
dme = 3*gaussFun(alpha, deltaT, t).';

numberOfRawSamples = floor(paramsSignal.Freqsamp*paramsSignal.Intetime);
delay = 10e-6;
totalSamples = numberOfRawSamples;
numberOfZeros = numberOfRawSamples - length(dme);

mixtureSignal = zeros(totalSamples, length(params.JNRVector), length(nbits), monteCarloLoops);
xHat = zeros(totalSamples, 2, length(params.JNRVector), length(nbits), monteCarloLoops);

for loopIndex = 1:monteCarloLoops
    loopIndex
    paramsSignal.numberOfGPSSignals = 1;
    paramsSignal.FreqDopp = 0;
    
    %------------------Using different GPS and noise signals
    
    Timeofthisloop = 0:totalSamples-1;
    Carrphase = mod(2*pi*(paramsSignal.FreqDopp)*Timeofthisloop/paramsSignal.Freqsamp,2*pi);
    Carrier = exp(1i*Carrphase).';
    
    interferenceSignal = [zeros(floor(numberOfZeros/4), 1);dme;zeros(floor(numberOfZeros*3/4), 1)];
    interferenceSignal = [interferenceSignal;zeros(length(Carrier) - length(interferenceSignal), 1)];
%     interferenceSignal = interferenceSignal.*Carrier;
    interferenceSignalPower = pow_eval(interferenceSignal);
    
    paramsSignal.FreqDopp = 0;
    GPSSignals = GPSGen(paramsSignal);
    GPSSignals = GPSSignals(1:numberOfRawSamples,:);
    GPSSignals = [GPSSignals(end - round(delay*params.fs)+1:end,:);GPSSignals(1:end - round(delay*params.fs),:)]; % Introducing artificial code delay
    GPSSignalsPower = pow_eval(GPSSignals);
    
%     noise = randn(length(GPSSignals), 1) + 1j*randn(length(GPSSignals), 1);
    noise = randn(length(GPSSignals), 1);
    noisePower = pow_eval(noise);

    for nbitsIndex = 1:length(nbits)
        nbitsIndex
        
        params.W0 = W0(:,:,nbitsIndex);
        for JNRIndex = 1:length(params.JNRVector)
            GPSSignalsAux = GPSSignals;
            interferenceSignalAux = interferenceSignal;
            GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
            mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
            interferenceSignalAux = interferenceSignalAux*sqrt(noisePower*10^(params.JNRVector(JNRIndex)/10)/interferenceSignalPower);
            mixtureSignal(:,JNRIndex,nbitsIndex,loopIndex) = mixtureGPS + interferenceSignalAux;
            mixtureSignal(:,JNRIndex,nbitsIndex,loopIndex) = quantise_gps(mixtureSignal(:,JNRIndex,nbitsIndex,loopIndex), nbits(nbitsIndex), 1);
        end
        
        [WTest, HTest, errorTest, PxxTest, f, t] = nmf_eval_v2(mixtureSignal(:,:,nbitsIndex,loopIndex), params);
        S = zeros([size(PxxTest{1,1}) 2]);
        
        for JNRIndex = 1:length(params.JNRVector)
            for i = 1:2
                S(:,:,i) = (W0(:,(i*params.numberOfSources/2) - (params.numberOfSources/2 -1):(i*params.numberOfSources/2),nbitsIndex) * ...
                    HTest{1,JNRIndex}((i*params.numberOfSources/2) - (params.numberOfSources/2 -1):(i*params.numberOfSources/2),:)./ (W0(:,:,nbitsIndex)*HTest{1,JNRIndex})).*PxxTest{1,JNRIndex};
                
                xHat(:,i,JNRIndex, nbitsIndex,loopIndex) = istft(S(:,:,i), params.fs, 'Window', params.window, 'OverlapLength', params.overlap, 'FFTLength', params.nfft);
                
            end
        end
    end
end

save(['.' filesep 'data' filesep 'nmf_testing_19.mat'], 'xHat', 'mixtureSignal', 'nbits', 'JNRVector');

rmpath(['.' filesep 'data']);
rmpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);