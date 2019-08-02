clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')


addpath(['..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep 'signalsGeneration' filesep]);

load sim_params_1.mat;

fs = paramsSignal.Freqsamp;

numberOfSources = 1;
random_state = 42;

params.fs = fs;
params.nfft = 128;
params.nperseg = 128;
params.overlap = params.nperseg-1;
params.hop_size = params.nperseg - params.overlap;
params.numberOfSources = numberOfSources;
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 10000;
params.tolChange = 1e-6;
params.tolError = 1e-6;
params.repetitions = 1;
SNR = -25;
params.JNRVector = [0];

numberOfRawSamples = 4096;
silenceSamples = round(20e-6*fs);

[interferenceSignal, GPSSignals] = signalGen(paramsSignal);
GPSSignals = GPSSignals(1:numberOfRawSamples,:);
interferenceSignal = interferenceSignal(1:numberOfRawSamples);

GPSSignals = [zeros(silenceSamples, size(GPSSignals, 2)); GPSSignals; zeros(silenceSamples, size(GPSSignals, 2))];
interferenceSignal = [zeros(silenceSamples, 1); interferenceSignal; zeros(silenceSamples, 1)];

interferenceSignalPower = pow_eval(interferenceSignal);
GPSSignalsPower = pow_eval(GPSSignals);

noise = randn(numberOfRawSamples + silenceSamples*2, 1) + 1j*randn(numberOfRawSamples + silenceSamples*2, 1);
noisePower = pow_eval(noise);

GPSMultiplier = sqrt(noisePower*10.^(SNR/10)./GPSSignalsPower);
mixtureGPS = sum(GPSSignals.*GPSMultiplier, 2) + noise;
interferenceSignal = interferenceSignal*sqrt(noisePower*10^(params.JNRVector/10)/interferenceSignalPower);

mixtureSignal = mixtureGPS + interferenceSignal;

if isreal(mixtureSignal)
    [PxxAux, f, t] = spectrogram(mixtureSignal, hann(params.nperseg), params.overlap, params.nfft, params.fs, 'power');
else
    [PxxAux, f, t] = spectrogram(mixtureSignal, hann(params.nperseg), params.overlap, params.nfft, params.fs, 'centered', 'power');
end

figure;
surf(t*1e6, f/1e6, 10*log10(abs(PxxAux).^2), 'EdgeColor', 'none');
axis xy;
view(0, 90);
xlim([min(t) max(t)]*1e6);
ylim([min(f) max(f)]/1e6);
ylabel('Frequency [MHz]');
xlabel('Time [$\mu$s]');
ax = gca;
set(ax, 'CLim', [-10 50], 'colormap', jet);
c = colorbar;
c.Label.String = '[dB]';

rmpath(['..' filesep '.' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep]);