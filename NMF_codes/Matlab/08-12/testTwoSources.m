clear;
clc;
close all;

addpath(['..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep '.' filesep 'Sigtools' filesep 'NMF_algorithms'])

fs = 32.768e6;
numberOfSources = 1;
secondsOfData = 8.62e-6;
secondsOfSilence = 100e-6;
numberOfSamples = secondsOfData*fs;
totalSamples = 4096;
bandwidth = 1e6;
f0 = 0;
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
params.tolChange = 1e-9;
params.tolError = 1e-9;
params.repetitions = 10;
params.JNRVector = [0];

rng(random_state);

%signal mixture definition---------
t = 0:1/fs:(secondsOfData - 1/fs);
f = ((bandwidth/2)/secondsOfData)*t + f0;
f1 = 1;

signal1 = exp(1j*2*pi*f1*f.*t).';
signal1 = repmat(signal1, ceil(100e-6/secondsOfData), 1);

f = ((bandwidth/2)/secondsOfData)*t + 5e6;

signal2 = exp(1j*2*pi*f1*f.*t).';
signal2 = repmat(signal2, ceil(100e-6/secondsOfData), 1);
onsetTime = 20e-6;
offsetTime = onsetTime + length(signal1)/fs;

signal1 = [zeros(round(onsetTime*fs), 1); signal1;zeros(round(onsetTime*fs), 1)];
signal2 = [zeros(round(onsetTime*fs/2), 1); signal2;zeros(round(onsetTime*fs*3/2), 1)];
signal1Length = length(signal1);

onset = find(signal1, 1, 'first') - params.nperseg;
offset = find(signal1, 1, 'last') - params.nperseg;

mixtureSignal = signal1 + signal2(1:end-1);
%--------------------------------------------

[W, ~,~, PxxAux, f, t] = nmf_eval(mixtureSignal, params);
plot(f/1e6, W{1, 1}(:,1))
figure;
findpeaks(W{1, 1}(:,1))

figure;
surf(t*1e6, f/1e6, 10*log10(abs(PxxAux{1,1}).^2), 'EdgeColor', 'none');
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

params.numberOfSources = 2;
[W, H, ~, PxxAux, f, t] = nmf_eval(mixtureSignal, params);

figure;
plot(f/1e6, W{1,1}(:,1))
figure;
plot(f/1e6, W{1,1}(:,2))

figure;
plot(t*1e6, H{1,1}(1,:))
figure;
plot(t*1e6, H{1,1}(2,:))

rmpath(['..' filesep '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep '.' filesep 'Sigtools' filesep])