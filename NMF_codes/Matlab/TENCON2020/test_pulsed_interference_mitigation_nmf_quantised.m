clear;
clc;
close all;

addpath(['..' filesep 'Sigtools' filesep])
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
addpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['.' filesep 'data']);

load nmf_training_quantised_01.mat;
load sim_params_3.mat;

monteCarloLoops = 1;
SNR = -25;

params.JNRVector = [10 20 30];
params.fs = paramsSignal.Freqsamp;
params.nfft = 256;
params.nperseg = 256;
params.overlap = params.nperseg - 1;
params.hop_size = params.nperseg - params.overlap;
params.window = kaiser(params.nperseg, 5);
params.specType = 'power';
params.numberOfSources = 50;
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 500;
params.tolChange = 1e-3;
params.tolError = 1e-3;
params.repetitions = 1;
params.verbose = true;
JNRVector = params.JNRVector;

numberOfRawSamples = floor(paramsSignal.Freqsamp*paramsSignal.Intetime);
delay = 10e-6;
totalSamples = numberOfRawSamples;

paramsSignal.numberOfGPSSignals = 1;
paramsSignal.FreqDopp = 3e3;
params.numberOfSources = params.numberOfSources*2;
params.init = 'custom';
params.W0 = W0;
params.transform = false;
params.alg = 'snmf';

paramsSignal.FreqDopp = 1e3;
GPSSignals = GPSGen(paramsSignal);
GPSSignals = GPSSignals(1:numberOfRawSamples,:);

GPSSignals = [GPSSignals(end - round(delay*params.fs)+1:end,:);GPSSignals(1:end - round(delay*params.fs),:)]; % Introducing artificial code delay
GPSSignalsPower = pow_eval(GPSSignals);

bandwidthVector = 30*1e6;
periodVector = 20e-6;
frontEndBandwidth = 2e6;

paramsSignal.Noneperiod = round(periodVector*params.fs);                   % number of samples with a sweep time
paramsSignal.Initphase = 0;
paramsSignal.IFmin = -bandwidthVector/2;                                                  % start frequency
paramsSignal.IFmax = bandwidthVector/2;                    % end frequency
paramsSignal.foneperiod(1:paramsSignal.Noneperiod) = linspace(paramsSignal.IFmin, paramsSignal.IFmax, paramsSignal.Noneperiod);
% paramsSignal.Intenumb = paramsSignal.Noneperiod;
[interferenceSignal, ~] = interferenceGen(paramsSignal);

interferenceSignal = [interferenceSignal;zeros(round(periodVector*params.fs), 1)];
interferenceSignal = repmat(interferenceSignal, 30, 1);
interferenceSignal = interferenceSignal(1:numberOfRawSamples);

Timeofthisloop = 0:length(interferenceSignal)-1;
Carrphase = mod(2*pi*(paramsSignal.FreqDopp)*Timeofthisloop/paramsSignal.Freqsamp,2*pi);
Carrier = exp(1i*Carrphase).';

interferenceSignal = interferenceSignal.*Carrier;
interferenceSignalPower = pow_eval(interferenceSignal);

[b,a] = butter(11, (frontEndBandwidth/(params.fs/2)));

nbitsVector = [2 4];
mixtureSignal = zeros(totalSamples, length(params.JNRVector), length(nbitsVector), 1, monteCarloLoops);
xHat = zeros(length(mixtureSignal), 2, length(params.JNRVector), length(nbitsVector), monteCarloLoops);

for loopIndex = 1:monteCarloLoops
    loopIndex
    
    %------------------Using different GPS and noise signals
    %     interferenceSignal = [zeros(floor(numberOfZeros/4), 1);dme;zeros(floor(numberOfZeros*3/4), 1)];
    %     interferenceSignal = [interferenceSignal;zeros(length(Carrier) - length(interferenceSignal), 1)];
    %     interferenceSignal = interferenceSignal.*Carrier;
    %     interferenceSignalPower = pow_eval(interferenceSignal);
    
    noise = randn(length(GPSSignals), 1) + 1j*randn(length(GPSSignals), 1);
    noisePower = pow_eval(noise);
    for nbitsIndex = 1:length(nbitsVector)
        params.W0 = W0{1,nbitsIndex};
        for JNRIndex = 1:length(params.JNRVector)
            GPSSignalsAux = GPSSignals;
            interferenceSignalAux = interferenceSignal;
            GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
            mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
            interferenceSignalAux = interferenceSignalAux*sqrt(noisePower*10^(params.JNRVector(JNRIndex)/10)/interferenceSignalPower);
            mixtureSignalAux = mixtureGPS + interferenceSignalAux;
            mixtureSignalAux = filter(b, a, mixtureSignalAux);

            mixtureSignal(:,JNRIndex,nbitsIndex,1,loopIndex) = quantise_gps(mixtureSignalAux, nbitsVector(nbitsIndex), 2);
            %         mixtureSignal(:,JNRIndex,1,1,loopIndex) = mixtureSignalAux;
        end
        
        [WTest, HTest, errorTest, PxxTest, f, t] = nmf_eval_v2(mixtureSignal(:,:,nbitsIndex,1,loopIndex), params);
        S = zeros([size(PxxTest{1,1}) 2]);
        
        for JNRIndex = 1:length(params.JNRVector)
            for i = 1:2
                S(:,:,i) = (params.W0(:,(i*params.numberOfSources/2) - (params.numberOfSources/2 -1):(i*params.numberOfSources/2)) * ...
                    HTest{1,JNRIndex}((i*params.numberOfSources/2) - (params.numberOfSources/2 -1):(i*params.numberOfSources/2),:)./ (params.W0*HTest{1,JNRIndex})).*PxxTest{1,JNRIndex};
                
                xHat(:,i,JNRIndex,nbitsIndex,loopIndex) = istft(S(:,:,i), params.fs, 'Window', params.window, 'OverlapLength', params.overlap, 'FFTLength', params.nfft);
                
            end
        end
    end
end

save(['.' filesep 'data' filesep 'nmf_testing_quantised_02.mat'], 'xHat', 'mixtureSignal', 'JNRVector', 'nbitsVector');

rmpath(['.' filesep 'data']);
rmpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);