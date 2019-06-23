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

signalLength = 125e-6;
numberOfSamples = round(signalLength*fs);
desiredSignalPower = db2pow(10);
monteCarloLoops = 500;


stftSignal = zeros(params.nfft, (numberOfSamples - params.nperseg + 1)/(params.nperseg - params.overlap), monteCarloLoops);

for i = 1:monteCarloLoops
    
    %signal mixture definition---------
    signal = randn(numberOfSamples, 1) + 1j*randn(numberOfSamples, 1);
    signalPower = signal'*signal/numberOfSamples;
    
    signal = signal*sqrt(desiredSignalPower/signalPower);
    mixtureSignal = signal;
    %--------------------------------------------
    
    mixtureSignalBuffered = buffer(mixtureSignal, params.nperseg, params.overlap, 'nodelay');
    window = diag(hann(params.nperseg));
    
    mixtureSignalBuffered = window*mixtureSignalBuffered;
    dftMatrix = dftmtx(params.nfft);
    
    stftSignal(:,:,i) = dftMatrix*mixtureSignalBuffered;
    stftSignal(:,:,i) = fftshift(stftSignal(:,:,i), 1);    
end

x = real(reshape(stftSignal, size(stftSignal,1)*size(stftSignal,3),size(stftSignal,2)));

figure
[N, edges] = histcounts(x(:,1), 100, 'Normalization', 'pdf');
edges = (edges(1:end-1) + edges(2:end))/2;
plot(edges, N);
hold on
pd_output = fitdist(x(:,1), 'Normal');
y = pdf(pd_output,edges);
plot(edges, y, '--');
legend('Data pdf', 'Fitted pdf');
xlabel('$s$');
ylabel('$f(s)$')

figure
[N, edges] = histcounts(x(:,1), 100, 'Normalization', 'cdf');
edges = (edges(1:end-1) + edges(2:end))/2;
plot(edges, N);
hold on
y = cdf(pd_output, edges);
plot(edges, y, '--');
xlabel('$s$');
ylabel('$F(s)$')
legend('Fitted', 'cdf');

[hreal,preal,kreal,creal] = lillietest(x(:,1)); 

rmpath(['..' filesep '.' filesep 'Sigtools' filesep])