clear;
clc;
close all;

addpath(['..' filesep 'Sigtools' filesep])
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
addpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])

load sim_params_1.mat;

monteCarloLoops = 100;
SNR = -25;
alpha = 4.5e11;
deltaT = 12e-6;

params.fs = paramsSignal.Freqsamp;
params.nfft = 256;
params.nperseg = 256;
params.overlap = params.nperseg - 1;
params.hop_size = params.nperseg - params.overlap;
params.window = ones(params.nperseg, 1);
params.specType = 'power';
params.numberOfSources = 1;
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 10000;
params.tolChange = 1e-6;
params.tolError = 1e-6;
params.repetitions = 1;
params.JNRVector = -15;

t = -10e-6 - 1/params.fs:1/params.fs:20e-6 - 1/params.fs;
gaussFun = @(alpha, deltaT, t) exp((-alpha/2).*t.^2) + exp((-alpha/2).*(t - deltaT).^2);
dme =  gaussFun(alpha, deltaT, t).';

numberOfRawSamples = 4096;
totalSamples = numberOfRawSamples;
numberOfZeros = numberOfRawSamples - length(dme);

interferenceSignal = [zeros(numberOfZeros/2, 1);dme;zeros(numberOfZeros/2, 1)];
Timeofthisloop = 0:totalSamples-1;               
Carrphase = mod(2*pi*(paramsSignal.FreqDopp)*Timeofthisloop/paramsSignal.Freqsamp,2*pi);
Carrier = exp(1i*Carrphase).';

interferenceSignal = interferenceSignal.*Carrier;
interferenceSignalPower = pow_eval(interferenceSignal);

GPSSignals = GPSGen(paramsSignal);
GPSSignals = GPSSignals(1:numberOfRawSamples,1);
GPSSignalsPower = pow_eval(GPSSignals);

for loopIndex = 1:monteCarloLoops
    loopIndex
    mixtureSignal = zeros(totalSamples, length(params.JNRVector));
    noise = randn(totalSamples, 1) + 1j*randn(totalSamples, 1);
    noisePower = pow_eval(noise);

    for i = 1:length(params.JNRVector)
        GPSSignalsAux = GPSSignals;
        interferenceSignalAux = interferenceSignal;
        GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
        mixtureGPS = sum(GPSSignalsAux.*GPSMultiplier, 2) + noise;
        interferenceSignalAux = interferenceSignalAux*sqrt(noisePower*10^(params.JNRVector(i)/10)/interferenceSignalPower);
        mixtureSignal(:,i) = mixtureGPS + interferenceSignalAux;
    end
    
    [W, ~, ~, Pxx, f, t] = nmf_eval_v2(mixtureSignal, params);
            
    for JNRIndex = 1:length(params.JNRVector)
        
        inputNMF = abs(Pxx{1, JNRIndex}).^2;
        figure;
        surf(t*1e6, f/1e6, 10*log10(inputNMF), 'EdgeColor', 'none');
        axis xy;
        view(0, 90);
        xlabel('Time [$\mu$s]');
        ylabel('Frequency [MHz]');
        xlim([0 max(t)]*1e6);
        ylim([min(f) max(f)]/1e6);
        ax = gca;
        set(ax, 'CLim', [0 40], 'colormap', jet);
        c = colorbar;
        c.Label.String = '[dB]';
        c.TickLabelInterpreter = 'latex';
        
        inputNMF = inputNMF - mean(inputNMF);
        inputNMF = inputNMF.*sqrt(1./var(inputNMF));
        inputNMFAux = sqrt(sum(inputNMF.*inputNMF)) + eps;
        inputNMFNormalised = inputNMF./inputNMFAux;
        
        WNormalised = W{1, JNRIndex}(:,1) - mean(W{1, JNRIndex}(:,1));
        WNormalised = WNormalised.*sqrt(1./var(WNormalised));
        WNormalised = WNormalised ./ (norm(WNormalised) + eps);
        
        output = inputNMFNormalised.'*WNormalised;
        figure;
        plot(W{1,1});
        figure;
        plot(output);
        
        for thresholdIndex = 1:length(thresholdVector)
            for window_median_length_index = 1:length(window_median_length_vector)
                detection_res(loopIndex, bandwidthIndex, periodIndex, JNRIndex, thresholdIndex, window_median_length_index) = ...
                    median(detection_eval(output, thresholdVector(thresholdIndex)));
            end
        end
    end

end

rmpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);