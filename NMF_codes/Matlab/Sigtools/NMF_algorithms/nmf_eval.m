
function [W, H, reconstructError, PxxAux, f, t] = nmf_eval(mixtureSignal)

% set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
% set(groot, 'defaultLegendInterpreter','latex');

addpath(['.' filesep 'Sigtools' filesep])
addpath(['.' filesep 'Sigtools' filesep 'NMF_algorithms'])


% figPath = './figs-02-14/';
fs = 32.768e6;
nfft = 256;

params.numberOfSources = 2;
params.init = 'random';
params.betaDivergence = 'frobenius';
params.numberOfIterations = 10000;
params.random_state = 42;
params.tolChange = 1e-6;
params.tolError = 1e-6;


% exp_name = 'TK';    
% f0 = 0;
% secondsOfData = 8.62e-6;
% numberOfSamples = secondsOfData*fs;
% totalSamples = 4096;
% bandwidth = 1e6;
JNRVector = [-10, -5, 0, 10];
JNRVector = inf;


% t = 0:1/fs:(secondsOfData - 1/fs);
% % f = ((bandwidth/2)/secondsOfData)*t + f0;
% f = ((bandwidth/2)/secondsOfData) + f0;
% f1 = 1;
% f2 = 4;
% 
% signal1 = sin(2*pi*f1*f.*t).';    
% signal1 = repmat(signal1, ceil(100e-6/secondsOfData), 1);
% signal2 = sin(2*pi*f2*f.*t).';  
% signal2 = repmat(signal2, ceil(100e-6/secondsOfData), 1);
% signal1Length = length(signal1);
% 
% signal1 = [zeros(totalSamples - signal1Length, 1); signal1];
% signal2 = [signal2;zeros(totalSamples - signal1Length, 1);];
% signal1Length = totalSamples;
% mixtureSignal = signal1 + signal2;
powMixture = pow_eval(mixtureSignal);
signal1Length = length(mixtureSignal);

noise = randn(signal1Length, 1) + 1j*randn(signal1Length, 1);
powNoise = pow_eval(noise);

for i = 1:length(JNRVector)
    
    powAux = powMixture/JNRVector(i);
    noise2 = noise*sqrt(powAux/powNoise);
    data = mixtureSignal + noise2;
    
    if all(isreal(data))
        [PxxAux, f, t] = spectrogram(data, hann(nfft), nfft-1, nfft, fs, 'power');
    else
        [PxxAux, f, t] = spectrogram(data, hann(nfft), nfft-1, nfft, fs, 'centered', 'power');
    end    
    
    dataCell{1} = PxxAux;
%     figure;
%     surf(t*1e6, f/1e6, 10*log10(abs(dataCell{1}).^2), 'EdgeColor', 'none');
%     axis xy;
%     axis tight;
%     colormap(jet); 
%     view(0,90);
%     ylabel('Frequency [MHz]');
%     xlabel('Time [\mu s]');
%     colorbar;
%     caxis([-120 20]);
    
%     dataCell{2} = TK_filtering(PxxAux);
    
    for j = 1:length(dataCell)
        inputNMF = abs(dataCell{j}).^2;
        [W, H, reconstructError] = nmf_py3(inputNMF, params);
        figure;
        plot(f, W);
        figure;
        plot(t, H);
        figure;
        plot(reconstructError(reconstructError~=0));
        
        for k = 1:params.numberOfSources
            separatedSignal = W(:,k)*H(k,:);
            figure;
            surf(t*1e6, f/1e6, 10*log10(abs(separatedSignal).^2), 'EdgeColor', 'none');
            axis xy;
            axis tight;
            colormap(jet);
            view(0,90);
            ylabel('Frequency [MHz]');
            xlabel('Time [\mu s]');
            colorbar;
        end
        
    end
        
end

rmpath(['.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['.' filesep 'Sigtools' filesep])