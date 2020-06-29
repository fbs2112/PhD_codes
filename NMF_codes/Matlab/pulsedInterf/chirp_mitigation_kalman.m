clear;
clc;
close all;

addpath(['..' filesep 'Sigtools' filesep])
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
addpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['.' filesep 'data']);
addpath(['.' filesep 'pai_fun']);

load sim_params_3.mat;

monteCarloLoops = 100;
SNR = -25;
nbits = 0;

params.fs = paramsSignal.Freqsamp;
params.JNRVector = [-5 0 10 30 50];
JNRVector = params.JNRVector;

numberOfRawSamples = floor(paramsSignal.Freqsamp*paramsSignal.Intetime);
delay = 10e-6;
totalSamples = numberOfRawSamples;

initialFrequency = 2e6;
bandwidthVector = (2:3:14)*1e6;
periodVector = 8.62e-6;

paramsSignal.Initphase = 0;
paramsSignal.Noneperiod = round(periodVector*params.fs);                   % number of samples with a sweep time

paramsSignal.FreqDopp = 1e3;

Timeofthisloop = 0:totalSamples-1;
Carrphase = mod(2*pi*(paramsSignal.FreqDopp)*Timeofthisloop/paramsSignal.Freqsamp,2*pi);
Carrier = exp(1i*Carrphase).';

paramsSignal.numberOfGPSSignals = 1;
GPSSignals = GPSGen(paramsSignal);
GPSSignals = GPSSignals(1:numberOfRawSamples,:);
GPSSignals = [GPSSignals(end - round(delay*params.fs)+1:end,:);GPSSignals(1:end - round(delay*params.fs),:)]; % Introducing artificial code delay
GPSSignalsPower = pow_eval(GPSSignals);

xHatPai = zeros(totalSamples, length(params.JNRVector), length(nbits), length(bandwidthVector), monteCarloLoops);
IFEstimationFlag = true;

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
    iflaw = zeros(totalSamples, length(bandwidthVector));
    
    for bandwidthIndex = 1:length(bandwidthVector)
        bandwidthIndex
        paramsSignal.IFmin = -bandwidthVector(bandwidthIndex)/2;                                                  % start frequency
        paramsSignal.IFmax = bandwidthVector(bandwidthIndex)/2;                    % end frequency
        paramsSignal.foneperiod(1:paramsSignal.Noneperiod) = linspace(paramsSignal.IFmin, paramsSignal.IFmax, paramsSignal.Noneperiod);
        
        [interferenceSignal, iflaw(:,bandwidthIndex)] = interferenceGen(paramsSignal);
        interferenceSignal = interferenceSignal(1:numberOfRawSamples);
        interferenceSignal = interferenceSignal.*Carrier;
        interferenceSignalPower = pow_eval(interferenceSignal);
        
        for nbitsIndex = 1:length(nbits)
            
            for JNRIndex = 1:length(params.JNRVector)
                GPSSignalsAux = GPSSignals;
                interferenceSignalAux = interferenceSignal;
                GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
                mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
                interferenceSignalAux = interferenceSignalAux*sqrt(noisePower*10^(params.JNRVector(JNRIndex)/10)/interferenceSignalPower);
                mixtureSignal = mixtureGPS + interferenceSignalAux;
                
                [iflawestiTF, IFrate, linearflag] = IFEstimator_pai(mixtureSignal, iflaw(:,bandwidthIndex), IFEstimationFlag);
                xHatPai(:,JNRIndex,nbitsIndex,bandwidthIndex,loopIndex) = KF_pai(mixtureSignal, IFEstimationFlag, iflawestiTF, IFrate, linearflag);
            end
            
        end
    end
end

save(['.' filesep 'data' filesep 'resultsPai04.mat'], 'xHatPai', 'nbits', 'JNRVector');

rmpath(['.' filesep 'pai_fun']);
rmpath(['.' filesep 'data']);
rmpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);