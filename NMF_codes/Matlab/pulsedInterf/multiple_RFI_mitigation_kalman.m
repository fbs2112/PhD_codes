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
SNR = -20;
nbits = 0;

params.fs = paramsSignal.Freqsamp;
params.JNRVector = (0:5:30);
JNRVector = params.JNRVector;

numberOfRawSamples = floor(paramsSignal.Freqsamp*paramsSignal.Intetime);
delay = 10e-6;
totalSamples = numberOfRawSamples;

bandwidthVector = (2:6:14)*1e6;
periodVector = 8.62e-6;
numberOfZerosVector = params.fs*(1e-3)*[0];
paramsSignal.Initphase = 0;
paramsSignal.Noneperiod = round(periodVector*params.fs);                   % number of samples with a sweep time

GPSSignals = GPSGen(paramsSignal);
GPSSignals = GPSSignals(1:numberOfRawSamples,:);
GPSSignals = [GPSSignals(end - round(delay*params.fs)+1:end,:);GPSSignals(1:end - round(delay*params.fs),:)]; % Introducing artificial code delay
GPSSignalsPower = pow_eval(GPSSignals);

CWFrequency = 2e6;
paramsSignal.IFmin = CWFrequency;                                          
paramsSignal.IFmax = CWFrequency;
paramsSignal.foneperiod(1:paramsSignal.Noneperiod) = linspace(paramsSignal.IFmin, paramsSignal.IFmax, paramsSignal.Noneperiod);
[interferenceSignalCW, ~] = interferenceGen(paramsSignal);

xHatPaiPerf = zeros(totalSamples, length(params.JNRVector), length(numberOfZerosVector), length(nbits), length(bandwidthVector), monteCarloLoops);
xHatPaiIFest = zeros(totalSamples, length(params.JNRVector), length(numberOfZerosVector), length(nbits), length(bandwidthVector), monteCarloLoops);

IFEstimationFlag = true; %perfect IF estimation
IFEstimationFlag2 = false; %IF estimation

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
        interferenceSignal = interferenceSignal + interferenceSignalCW;
        
        for nbitsIndex = 1:length(nbits)
            for numberOfZerosIndex = 1:length(numberOfZerosVector)  
                for JNRIndex = 1:length(params.JNRVector)
                    GPSSignalsAux = GPSSignals;
                    interferenceSignalAux = interferenceSignal;
                    interferenceSignalAux(1:numberOfZerosVector(numberOfZerosIndex)) = 0;
                    interferenceSignalPower = pow_eval(interferenceSignalAux);
                    GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
                    mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
                    interferenceSignalAux = interferenceSignalAux*sqrt(noisePower*10^(params.JNRVector(JNRIndex)/10)/interferenceSignalPower);
                    mixtureSignal = mixtureGPS + interferenceSignalAux;

                    [iflawestiTF, IFrate, linearflag] = IFEstimator_pai(mixtureSignal, iflaw(:,bandwidthIndex), IFEstimationFlag);
                    xHatPaiPerf(:,JNRIndex,numberOfZerosIndex,nbitsIndex,bandwidthIndex,loopIndex) = KF_pai(mixtureSignal, IFEstimationFlag, iflawestiTF, IFrate, linearflag);

                    [iflawestiTF, IFrate, linearflag] = IFEstimator_pai(mixtureSignal, iflaw(:,bandwidthIndex), IFEstimationFlag2);
                    xHatPaiIFest(:,JNRIndex,numberOfZerosIndex,nbitsIndex,bandwidthIndex,loopIndex) = KF_pai(mixtureSignal, IFEstimationFlag2, iflawestiTF, IFrate, linearflag);
                end
            end
            
        end
    end
end

xHatPaiPerf = single(xHatPaiPerf);
xHatPaiIFest = single(xHatPaiIFest);

save(['.' filesep 'data' filesep 'resultsPai11.mat'], 'xHatPaiPerf', 'xHatPaiIFest', 'nbits', 'JNRVector', '-v7.3');

rmpath(['.' filesep 'pai_fun']);
rmpath(['.' filesep 'data']);
rmpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);