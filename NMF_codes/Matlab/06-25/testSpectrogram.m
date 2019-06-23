clear;
clc;
close all;

addpath(['..' filesep '.' filesep 'Sigtools' filesep])

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

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
params.tolChange = 1e-6;
params.tolError = 1e-6;
params.repetitions = 1;
params.JNRVector = [10];


%signal mixture definition---------
t = 0:1/fs:(secondsOfData - 1/fs);
f = ((bandwidth/2)/secondsOfData)*t + f0;
f1 = 1;

signal1 = exp(1j*2*pi*f1*f.*t).';
signal1 = repmat(signal1, ceil(100e-6/secondsOfData), 1);
onsetTime = 20e-6;
offsetTime = onsetTime + length(signal1)/fs;

signal1 = [zeros(round(onsetTime*fs), 1); signal1;zeros(round(onsetTime*fs), 1)];
mixtureSignal = signal1;
%--------------------------------------------

powMixture = pow_eval(mixtureSignal);
signal1Length = length(mixtureSignal);
noise = randn(signal1Length, 1) + 1j*randn(signal1Length, 1);
powNoise = pow_eval(noise);
powAux = powMixture/db2pow(params.JNRVector);
noise2 = noise*sqrt(powAux/powNoise);
mixtureSignal = mixtureSignal + noise2;

mixtureSignalBuffered = buffer(mixtureSignal, params.nperseg, params.overlap, 'nodelay');
window = diag(hann(params.nperseg));

mixtureSignalBuffered = window*mixtureSignalBuffered;
dftMatrix = dftmtx(params.nfft);

stftSignal = dftMatrix*mixtureSignalBuffered;
stftSignal = fftshift(stftSignal, 1);
inputNMF = abs(stftSignal).^2;
figure;
surf(10*log10(inputNMF), 'EdgeColor', 'none');
axis xy;
view(0, 90);
% xlim([min(t) max(t)]*1e6);
% ylim([min(f) max(f)]/1e6);
ylabel('Frequency [MHz]');
xlabel('Time [$\mu$s]');
ax = gca;
set(ax, 'CLim', [-60 20], 'colormap', jet);
c = colorbar;
c.Label.String = '[dB]';
c.TickLabelInterpreter = 'latex';


[hreal,preal,kreal,creal] = lillietest(real(stftSignal(:,1))); 
[himag,pimag,kimag,cimag] = lillietest(imag(stftSignal(:,1))); 
[hnoise,pnoise,knoise,cnoise] = lillietest(real(noise2)); 

Wvector = ones(params.nfft, 1);

[muhat, muci] = expfit(inputNMF(:,1));

figure;
data = linspace(min(inputNMF(:,1)), max(inputNMF(:,1)), 100);
plot(exppdf(data, muhat));

[htest,ptest,ktest,ctest] = lillietest(inputNMF(:,1), 'Distribution', 'exponential'); 



figure
[N, edges] = histcounts(x, 100, 'Normalization', 'pdf');
edges = (edges(1:end-1) + edges(2:end))/2;
plot(edges, N);
hold on
pd_output = fitdist(x, 'Normal');
y = pdf(pd_output,edges);
plot(edges, y, '--');
legend('Data pdf', 'Fitted pdf');
xlabel('$s$');
ylabel('$f(s)$')
xlim([-4 4]);

rmpath(['..' filesep '.' filesep 'Sigtools' filesep])