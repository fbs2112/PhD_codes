clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep '.' filesep 'Sigtools' filesep 'NMF_algorithms'])

fs = 32.768e6;
numberOfSources = 1;
secondsOfData = 8.62e-6;
secondsOfSilence = 100e-6;
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
params.JNRVector = [inf];

rng(random_state);


signalLength = 125e-6;
numberOfSamples = round(signalLength*fs);
desiredSignalPower = db2pow(10);

window_length = round(3e-6*fs);
window_median_length = 201; 
similarityName = 'inner';
stdVector = 0;

if strcmp(similarityName, 'gaussian')
    stdVector = 7:2:13;
end

monteCarloLoops = 200;

output = zeros((numberOfSamples - params.nperseg + 1)/(params.nperseg - params.overlap), length(params.JNRVector), monteCarloLoops);
% outputZscore = zeros((numberOfSamples - params.nperseg + 1)/(params.nperseg - params.overlap), length(params.JNRVector), monteCarloLoops);
% outputNotNorm = zeros((numberOfSamples - params.nperseg + 1)/(params.nperseg - params.overlap), length(params.JNRVector), monteCarloLoops);
% outputNormNotZscore = zeros((numberOfSamples - params.nperseg + 1)/(params.nperseg - params.overlap), length(params.JNRVector), monteCarloLoops);

WTensor = zeros(params.nfft, length(params.JNRVector), monteCarloLoops);

for loopIndex = 1:monteCarloLoops
    loopIndex
    
    %signal mixture definition---------
    signal = randn(numberOfSamples, 1) + 1j*randn(numberOfSamples, 1);
    signalPower = signal'*signal/numberOfSamples;
    
    signal = signal*sqrt(desiredSignalPower/signalPower);
    mixtureSignal = signal;
    %--------------------------------------------
    [W, H, reconstructError, PxxAux, f, t] = nmf_eval(mixtureSignal, params, false);
    
    for JNRIndex = 1:length(params.JNRVector)
        inputNMF = abs(PxxAux{1, JNRIndex}).^2;
        for stdIndex = 1:length(stdVector)
            WTensor(:, JNRIndex, loopIndex) = W{1, JNRIndex}(:,1);
            for tIndex = 1:length(t)
                if strcmp(similarityName, 'gaussian')
                    output(tIndex, JNRIndex, loopIndex) = similarity_eval(inputNMF(:,tIndex), W{1, JNRIndex}(:,1), similarityName, stdVector(stdIndex), true);
                else
                    output(tIndex, JNRIndex, loopIndex) = similarity_eval(inputNMF(:,tIndex), W{1, JNRIndex}(:,1), similarityName, 'normalized', true);
%                     outputZscore(tIndex, JNRIndex, loopIndex) = similarity_eval(inputNMF(:,tIndex), W{1, JNRIndex}(:,1), similarityName, [], true);
%                     outputNotNorm(tIndex, JNRIndex, loopIndex) = similarity_eval(inputNMF(:,tIndex), W{1, JNRIndex}(:,1), similarityName, [], false);
%                     outputNormNotZscore(tIndex, JNRIndex, loopIndex) = similarity_eval(inputNMF(:,tIndex), W{1, JNRIndex}(:,1), similarityName, 'normalized', false);
                end
            end            
        end
    end
end

figure
histogram(WTensor(:,1,:), 50);
x = output(:,1,:);
x = x(:);

x = x - mean(x);
x = x*sqrt(1/var(x));
figure
[N, edges] = histcounts(x, 100, 'Normalization', 'pdf');
edges = (edges(1:end-1) + edges(2:end))/2;
plot(edges, N);
hold on
pd_output = fitdist(x, 'Normal');
y = pdf(pd_output,edges);
plot(edges, y);

figure
[N, edges] = histcounts(x, 100, 'Normalization', 'cdf');
edges = (edges(1:end-1) + edges(2:end))/2;
plot(edges, N);
hold on
y = cdf(pd_output, edges);
plot(edges, y);
[y, edge] = ecdf(x);
plot(edge, y);

legend('Fitted', 'cdf', ' ecdf');

[hjb,pjb,kjb,cjb] = jbtest(x); 
[hlillie,plillie,klillie,clillie] = lillietest(x); 

gaussian_fun = @(data, x, y) (1./sqrt(2*pi*y)) .* exp(-(data - x).^2/(2*y));
data = -10:0.1:10;
mean_gauss = 0;
var_gauss = 1;
gauss_data = gaussian_fun(data, mean_gauss, var_gauss);

figure
plot(data, gauss_data);
% [h,p,k,c] = jbtest(x);  

rmpath(['..' filesep '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep '.' filesep 'Sigtools' filesep])