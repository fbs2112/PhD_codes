clear;
clc;
close all;

addpath(['..' filesep '..' filesep 'Sigtools' filesep])
addpath(['..' filesep '..' filesep 'signalsGeneration' filesep]);

load(['..' filesep '..' filesep 'signalsGeneration' filesep 'sim_params' filesep 'sim_params_3.mat']);

monteCarloLoops = 1000;
SNR = -20;
numberOfSources = 2;

numberOfRawSamples = floor(paramsSignal.Freqsamp*paramsSignal.Intetime);
delay = 10e-6;
totalSamples = numberOfRawSamples;
paramsSignal.Initphase = 0;
paramsSignal.FreqDopp = 1e3;

GPSSignals = GPSGen(paramsSignal);
GPSSignals = GPSSignals(1:numberOfRawSamples,:);
GPSSignalsPower = pow_eval(GPSSignals);

for loopIndex = 1:monteCarloLoops
    loopIndex
    noise = randn(totalSamples, 1) + 1j*randn(totalSamples, 1);
    noisePower = pow_eval(noise);
    
    GPSSignalsAux = GPSSignals;
    
    GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
    mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
    mixtureSignal = mixtureGPS;
    wav_output(:,:,loopIndex) = wav_eval(mixtureSignal);
    
end
wav_output = abs(wav_output);
wav_output_permute = permute(wav_output,[1 3 2]);
wav_output_reshaped = reshape(wav_output,size(wav_output, 1)*size(wav_output, 3), size(wav_output, 2));

for i = 1:16
    figure;
    histogram(wav_output_reshaped(:,i), 20);
    varScales(i) = var(wav_output_reshaped(:,i));
end

varNoise = 1;
pfa = 0.01;
% PB_threshold = sqrt(-2*varNoise*log(pfa));


wav_output_reshaped_zeroed = wav_output_reshaped;
% wav_output_reshaped_zeroed(abs(wav_output_reshaped_zeroed) > PB_threshold) = 0;

for i = 1:16
   PB_threshold(i) = sqrt(-2*varNoise*log(pfa));
   aux = wav_output_reshaped_zeroed(:,i);
   aux(abs(aux) > PB_threshold(i)) = 0;
   zeroed_samples(i) = length(find(aux == 0))/length(aux);
end

rmpath(['..' filesep '..' filesep 'Sigtools' filesep])
rmpath(['..' filesep '..' filesep 'signalsGeneration' filesep]);