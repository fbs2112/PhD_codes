clear;
clc;
close all;

addpath(['..' filesep '..' filesep 'Sigtools' filesep])
addpath(['..' filesep '..' filesep 'signalsGeneration' filesep]);

load(['.' filesep '..' filesep 'data' filesep 'nmf_training_5.mat']);
load(['..' filesep '..' filesep 'signalsGeneration' filesep 'sim_params' filesep 'sim_params_3.mat']);

monteCarloLoops = 500;
SNR = -20;

JNRVector = 0:5:30;
varNoise = 1;
pfa = 0.01;
PB_threshold = sqrt(-2*varNoise*log(pfa));
numberOfZerosVector = params.fs*(1e-3)*[0];

numberOfRawSamples = floor(paramsSignal.Freqsamp*paramsSignal.Intetime);
delay = 10e-6;
totalSamples = numberOfRawSamples;

bandwidthVector = [2 8 14]*1e6;
periodVector = 8.62e-6;

paramsSignal.Noneperiod = round(periodVector*params.fs);                   % number of samples with a sweep time
paramsSignal.Initphase = 0;
paramsSignal.FreqDopp = 1e3;

GPSSignals = GPSGen(paramsSignal);
GPSSignals = GPSSignals(1:numberOfRawSamples,:);
GPSSignals = [GPSSignals(end - round(delay*params.fs)+1:end,:);GPSSignals(1:end - round(delay*params.fs),:)]; % Introducing artificial code delay
GPSSignalsPower = pow_eval(GPSSignals);

mixtureSignal = zeros(totalSamples, length(JNRVector), length(numberOfZerosVector), length(bandwidthVector), monteCarloLoops);
xHat = zeros(totalSamples, length(JNRVector), length(numberOfZerosVector), length(bandwidthVector), monteCarloLoops);

for loopIndex = 1:monteCarloLoops
    loopIndex
    
    if paramsSignal.FreqDopp ~= 0
        noise = randn(totalSamples, 1) + 1j*randn(totalSamples, 1);
        noiseVar = 2;
    else
        noise = randn(totalSamples, 1);
        noiseVar = 1;
    end
    noisePower = pow_eval(noise);
    
    for bandwidthIndex = 1:length(bandwidthVector)
        bandwidthIndex
        paramsSignal.IFmin = -bandwidthVector(bandwidthIndex)/2;                                                  % start frequency
        paramsSignal.IFmax = bandwidthVector(bandwidthIndex)/2;                    % end frequency
        paramsSignal.foneperiod(1:paramsSignal.Noneperiod) = linspace(paramsSignal.IFmin, paramsSignal.IFmax, paramsSignal.Noneperiod);
        
        [interferenceSignal, ~] = interferenceGen(paramsSignal);
        interferenceSignal = interferenceSignal(1:numberOfRawSamples);
        interferenceSignalPower = pow_eval(interferenceSignal);
        
        for numberOfZerosIndex = 1:length(numberOfZerosVector)
            for JNRIndex = 1:length(JNRVector)
                GPSSignalsAux = GPSSignals;
                interferenceSignalAux = interferenceSignal;
                interferenceSignalAux(1:numberOfZerosVector(numberOfZerosIndex)) = 0;
                GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
                mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
                interferenceSignalAux = interferenceSignalAux*sqrt(noisePower*10^(JNRVector(JNRIndex)/10)/interferenceSignalPower);
                mixtureSignal(:,JNRIndex,numberOfZerosIndex,bandwidthIndex,loopIndex) = mixtureGPS + interferenceSignalAux;
                
                wavOutput = wav_dec_eval(mixtureSignal(:,JNRIndex,numberOfZerosIndex,bandwidthIndex,loopIndex)); %wavelet decomposition
                wavOutputZeroed = wavOutput;
                wavOutputZeroed(abs(wavOutputZeroed) > PB_threshold) = 0;      %RFI mitigation 
                xHatAux = wav_rec_eval(wavOutputZeroed);                       %signal reconstruction
                xHat(:,JNRIndex,numberOfZerosIndex,bandwidthIndex,loopIndex) = xHatAux(1:numberOfRawSamples);
            end
        end
        
    end
end

xHat = single(xHat);

save(['..' filesep 'data' filesep 'wav_testing_1.mat'], 'xHat', 'mixtureSignal', 'JNRVector', '-v7.3');

rmpath(['..' filesep '..' filesep 'Sigtools' filesep])
rmpath(['..' filesep '..' filesep 'signalsGeneration' filesep]);