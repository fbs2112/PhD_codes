clear;
clc;
close all;

addpath(['..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep '.' filesep 'Misc'])


set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')


linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
dataPath = ['..' filesep '.' filesep 'figs' filesep '06-25' filesep];

fs = 32.768e6;
random_state = 42;
rng(random_state);
params.fs = fs;
params.nfft = 128;
params.nperseg = 128;
params.overlap = params.nperseg-1;
params.hop_size = params.nperseg - params.overlap;

signalLength = 125e-6;
numberOfSamples = round(signalLength*fs);
desiredSignalPower = db2pow(10);
monteCarloLoops = 100;

stftSignal = zeros(params.nfft, (numberOfSamples - params.nperseg + 1)/(params.nperseg - params.overlap), monteCarloLoops);
wVector = rand([params.nfft, 1]);

for i = 1:monteCarloLoops
    
    %signal mixture definition---------
    signal = randn(numberOfSamples, 1) + 1j*randn(numberOfSamples, 1);
    signalPower = signal'*signal/numberOfSamples;
    
    signal = signal*sqrt(desiredSignalPower/signalPower);
    mixtureSignal = signal;
    %--------------------------------------------
    
    mixtureSignalBuffered = buffer(mixtureSignal, params.nperseg, params.overlap, 'nodelay');
%     window = diag(hann(params.nperseg));
    window = eye(params.nperseg);
    
    mixtureSignalBuffered = window*mixtureSignalBuffered;
    dftMatrix = dftmtx(params.nfft);
    
    stftSignal(:,:,i) = dftMatrix*mixtureSignalBuffered;
    stftSignal(:,:,i) = fftshift(stftSignal(:,:,i), 1);    
    
    inputNMF2(:,:,i) = abs(stftSignal(:,:,i)).^2;
    s(:,i) = inputNMF2(:,:,i).'*wVector;
end

inputNMF = abs(stftSignal).^2;
% b = squeeze(sum(inputNMF, 1)).*wVector;

xNMF = reshape(inputNMF, size(inputNMF,1)*size(inputNMF,3),size(inputNMF,2));
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
xlabel('$\overline{y}$');
ylabel('$f(\overline{y})$')
axis tight;

figure
[N, edges] = histcounts(x(:,1), 100, 'Normalization', 'cdf');
edges = (edges(1:end-1) + edges(2:end))/2;
plot(edges, N);
hold on
y = cdf(pd_output, edges);
plot(edges, y, '--');
xlabel('$\overline{y}$');
ylabel('$F(\overline{y})$')
legend('Fitted', 'cdf');
axis tight;
[hreal,preal,kreal,creal] = lillietest(x(:,1)); 


figure
[N, edges] = histcounts(xNMF(:,1), 100, 'Normalization', 'pdf');
edges = (edges(1:end-1) + edges(2:end))/2;
plot(edges, N);
hold on
pd_output = fitdist(xNMF(:,1), 'exponential');
y = pdf(pd_output,edges);
plot(edges, y, '--');
legend('Data pdf', 'Fitted pdf');
xlabel('$x$');
ylabel('$f(x)$')
axis tight;

figure
[N, edges] = histcounts(xNMF(:,1), 100, 'Normalization', 'cdf');
edges = (edges(1:end-1) + edges(2:end))/2;
plot(edges, N);
hold on
y = cdf(pd_output, edges);
plot(edges, y, '--');
xlabel('$x$');
ylabel('$F(x)$')
legend('Fitted', 'cdf');
axis tight;

[hexpo,pexpo,kexpo,cexpo] = lillietest(xNMF(:,1), 'Distribution', 'exponential'); 

% wVector = rand([size(xNMF, 1), 1]);
% 
% for i = 1:size(xNMF, 2)
%     sVector(i) = wVector.'*xNMF(:,i);
% end

% b = sum(xNMF, 1);

pd = fitdist(xNMF(:,1), 'exponential');
hexp = chi2gof(xNMF(:,1),'CDF',pd);

for j = 1:size(s, 1)
    pd = fitdist(s(j,:).', 'Gamma');
    [hgamma(j), pgamma(j)] = chi2gof(s(j,:).', 'CDF', pd);

    pd = fitdist(s(j,:).', 'normal');
    [hgaussian(j), pgaussian(j)] = chi2gof(s(j,:).', 'CDF', pd);
end

figure;
stem(hgamma);

rateGamma = length(find(~hgamma))/length(hgamma);
rateGaussian = length(find(~hgaussian))/length(hgaussian);
figure;
stem(hgaussian);


figure
[N, edges] = histcounts(sVector.', 100, 'Normalization', 'pdf');
edges = (edges(1:end-1) + edges(2:end))/2;
plot(edges, N);
hold on
pd_output = fitdist(sVector.', 'Gamma');
y = pdf(pd_output,edges);
plot(edges, y, '--');
legend('Data pdf', 'Fitted pdf');
xlabel('$s$');
ylabel('$f(s)$')
axis tight;

figure
[N, edges] = histcounts(sVector.', 100, 'Normalization', 'cdf');
edges = (edges(1:end-1) + edges(2:end))/2;
plot(edges, N);
hold on
y = cdf(pd_output, edges);
plot(edges, y, '--');
xlabel('$s$');
ylabel('$F(s)$')
legend('Fitted', 'cdf');
axis tight;

rmpath(['..' filesep '.' filesep 'Misc'])
rmpath(['..' filesep '.' filesep 'Sigtools' filesep])