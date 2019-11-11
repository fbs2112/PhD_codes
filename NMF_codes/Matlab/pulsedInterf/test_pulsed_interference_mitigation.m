clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');
    
addpath(['..' filesep 'Sigtools' filesep])
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
addpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['..' filesep 'Misc'])

load sim_params_1.mat;

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
dataPath = ['..' filesep 'figs' filesep '11-12' filesep];

monteCarloLoops = 1;
SNR = -25;
alpha = 4.5e11;
deltaT = 12e-6;
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
params.numberOfIterations = 1000;
params.tolChange = 1e-6;
params.tolError = 1e-6;
params.repetitions = 1;

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
    
    plot(abs(interferenceSignal));
    ylabel('Magnitude');
    xlabel('Samples');
    axis tight;
    formatFig(gcf, [dataPath  'interf_time'], 'en', figProp);
    
    plot(abs(mixtureGPS));
    ylabel('Magnitude');
    xlabel('Samples');
    axis tight;
    formatFig(gcf, [dataPath  'gps_time'], 'en', figProp);

    [WInterf, HInterf, ~, PxxInterf, ~, ~] = nmf_eval_v2(interferenceSignal, params);
    
    [WSignal, HSignal, ~, PxxSignal, ~, ~] = nmf_eval_v2(mixtureGPS, params);
    
    W0 = [WInterf{1,1} WSignal{1,1}];
    params.numberOfSources = params.numberOfSources*2;
    params.init = 'custom';
    params.W0 = W0;
    params.transform = false;
    [~, H, ~, Pxx, f, t] = nmf_eval_v2(mixtureSignal, params);
    
    figure;
    surf(t*1e6, f/1e6, 10*log10(abs(PxxInterf{1,1} + 1e-6).^2), 'EdgeColor', 'none');
    axis xy;
    view(0, 90);
    xlabel('Time [$\mu$s]');
    ylabel('Frequency [MHz]');
    xlim([0 max(t)]*1e6);
    ylim([min(f) max(f)]/1e6);
    ax = gca;
    set(ax, 'CLim', [-20 40], 'colormap', jet);
    c = colorbar;
    c.Label.String = '[dB]';
    c.TickLabelInterpreter = 'latex';
    formatFig(gcf, [dataPath  'interf_spec'], 'en', figProp);
    
    figure;
    surf(t*1e6, f/1e6, 10*log10(abs(PxxSignal{1,1} + 1e-6).^2), 'EdgeColor', 'none');
    axis xy;
    view(0, 90);
    xlabel('Time [$\mu$s]');
    ylabel('Frequency [MHz]');
    xlim([0 max(t)]*1e6);
    ylim([min(f) max(f)]/1e6);
    ax = gca;
    set(ax, 'CLim', [-20 40], 'colormap', jet);
    c = colorbar;
    c.Label.String = '[dB]';
    c.TickLabelInterpreter = 'latex';
    formatFig(gcf, [dataPath  'gps_spec'], 'en', figProp);
    
    for i = 1:2
        S(:,:,i) = (W0(:,(i*params.numberOfSources/2) - (params.numberOfSources/2 -1):(i*params.numberOfSources/2)) * ...
            H{1,1}((i*params.numberOfSources/2) - (params.numberOfSources/2 -1):(i*params.numberOfSources/2),:)./ (W0*H{1,1})).*Pxx{1,1};
        
        xHat(:,i) = istft(S(:,:,i), params.fs, 'Window', ones(params.nperseg, 1), 'OverlapLength', params.overlap, 'FFTLength', params.nfft);
        
        figure;
        surf(t*1e6, f/1e6, 10*log10(abs(S(:,:,i)+ 1e-6).^2), 'EdgeColor', 'none');
        axis xy;
        view(0, 90);
        xlabel('Time [$\mu$s]');
        ylabel('Frequency [MHz]');
        xlim([0 max(t)]*1e6);
        ylim([min(f) max(f)]/1e6);
        ax = gca;
        set(ax, 'CLim', [-20 40], 'colormap', jet);
        c = colorbar;
        c.Label.String = '[dB]';
        c.TickLabelInterpreter = 'latex';
        formatFig(gcf, [dataPath  'spec_nmf_' num2str(i)], 'en', figProp);
    end
    
    figure;
    plot(abs(xHat(:,1)));
    ylabel('Magnitude');
    xlabel('Samples');
    axis tight;
    formatFig(gcf, [dataPath  'interf_time_nmf'], 'en', figProp);
    
    figure;
    plot(abs(xHat(:,2)));
    ylabel('Magnitude');
    xlabel('Samples');
    axis tight;
    formatFig(gcf, [dataPath  'gps_time_nmf'], 'en', figProp);
    
    [a,b] = xcorr(xHat(:,2), mixtureGPS);
    figure;
    plot(b, abs(a));
    ylabel('Magnitude');
    xlabel('Lags');
    axis tight;
    formatFig(gcf, [dataPath  'corr'], 'en', figProp);

end

rmpath(['..' filesep 'Misc'])
rmpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);