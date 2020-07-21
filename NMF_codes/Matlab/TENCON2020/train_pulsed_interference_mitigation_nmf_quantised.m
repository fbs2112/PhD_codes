clear;
clc;
close all;
    
addpath(['..' filesep 'Sigtools' filesep])
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
addpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])

load sim_params_3.mat;

monteCarloLoops = 1;
SNR = -25;
params.JNRVector = 0;

params.fs = paramsSignal.Freqsamp;
params.nfft = 256;
params.nperseg = 256;
params.overlap = params.nperseg - 1;
params.hop_size = params.nperseg - params.overlap;
% params.window = ones(params.nperseg, 1);
params.window = kaiser(params.nperseg, 2);

params.specType = 'power';
params.numberOfSources = 50;
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 500;
params.tolChange = 1e-3;
params.tolError = 1e-3;
params.repetitions = 1;
params.alg = 'snmf';
params.verbose = true;
% paramsSignal.Intenumb = round(100e-6*params.fs);
numberOfRawSamples = paramsSignal.Intenumb;
totalSamples = paramsSignal.Intenumb;

paramsSignal.numberOfGPSSignals = 1; %overriding number of GPS signals from configuration file
GPSSignals = GPSGen(paramsSignal);
GPSSignals = GPSSignals(1:numberOfRawSamples,:);
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

% interferenceSignal = [interferenceSignal;zeros(round(periodVector*params.fs), 1)];
interferenceSignal = repmat(interferenceSignal, 30, 1);
interferenceSignal = interferenceSignal(1:numberOfRawSamples);

Timeofthisloop = 0:length(interferenceSignal)-1;               
Carrphase = mod(2*pi*(paramsSignal.FreqDopp)*Timeofthisloop/paramsSignal.Freqsamp,2*pi);
Carrier = exp(1i*Carrphase).';

interferenceSignal = interferenceSignal.*Carrier;
interferenceSignalPower = pow_eval(interferenceSignal);
nbitsVector = [2 4];
[b,a] = butter(11, (frontEndBandwidth/(params.fs/2)));
% [h,w] = freqz(b,a,2^10,'whole',params.fs);
% interferenceSignal = interferenceSignal(1:round(150e-6*params.fs));

% interferenceSignal = buffer(interferenceSignal,params.nperseg, 0, 'nodelay');
% interferenceSignal = filter(b, a, interferenceSignal);
% interferenceSignal = interferenceSignal(:);
% % interferenceSignal = interferenceSignal(1:round(60e-6*params.fs));
% plot(w,(20*log10(abs(h))));
for loopIndex = 1:monteCarloLoops
    W0 = cell(1,length(nbitsVector));
    for nbitsIndex = 1:length(nbitsVector)
        loopIndex
        noise = randn(totalSamples, 1) + 1j*randn(totalSamples, 1);
        noisePower = pow_eval(noise);

        GPSSignalsAux = GPSSignals;
        GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
        mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
        mixtureGPS = filter(b, a, mixtureGPS);
        mixtureGPS = quantise_gps(mixtureGPS, nbitsVector(nbitsIndex), 2);
        
        interferenceAux = filter(b, a, interferenceSignal*3);
        interferenceAux = quantise_gps(interferenceAux, nbitsVector(nbitsIndex), 2);

        [WInterf, HInterf, errorInterfTrain, PxxInterf, ~, ~] = nmf_eval_v2(interferenceAux, params);
        [WSignal, HSignal, errorSignalTrain, PxxSignal, ~, ~] = nmf_eval_v2(mixtureGPS, params);

        W0{1,nbitsIndex} = [WInterf{1,1} WSignal{1,1}];
    end
end

save(['.' filesep 'data' filesep 'nmf_training_quantised_01.mat'], 'W0');

rmpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);