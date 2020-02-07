clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');

addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
addpath(['..' filesep 'Misc'])
addpath(['..' filesep 'Sigtools' filesep])

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
dataPath = ['..' filesep 'figs' filesep '02-06-2020' filesep];

load sim_params_2.mat;

load(['.' filesep 'data' filesep 'nmf_testing_12.mat']);

varN = 2;
Pfa = 0.05;
thresholdBorio = sqrt(-2*varN*log(Pfa));
fi = 0;
fs = paramsSignal.Freqsamp;
Nd = 81;
DopStep = 125; %(Hz)

trueDelay = 10e-6;
trueDoppler = 1.5e3;
paramsSignal.numberOfGPSSignals = 1;

[~, locC] = GPSGen(paramsSignal);   % Resample the code at data sampling frequency

JNRVector = [-5 10 30 50];
JNRVector = [0];


Tc = paramsSignal.Intetime;     % Coherent integration time in ms
Nc = floor(Tc * fs);

K = 1;  % Number of non-coherent integrations
codeDl = (0:(Nc - 1)) / fs; % Code delays

nbits = [2 4 8 16];
nbits = 0;
% for nbitsIndex = 1:length(nbits)
for nbitsIndex = 1:1
%     locC = quantise_gps(locC, nbits(nbitsIndex));
    for JNRIndex = 1:length(JNRVector)
        
        xHatBorio = mixtureSignal(:,JNRIndex,nbitsIndex);
        zeroedSamples(JNRIndex) = length(find(abs(xHatBorio) > thresholdBorio));
        fprintf('Number of zeroed samples: %i\n', zeroedSamples(JNRIndex));
        xHatBorio(abs(xHatBorio) > thresholdBorio) = 0;
        
        GPSSignals(1,:) = xHat(:,2,JNRIndex,nbitsIndex).';
        GPSSignals(2,:) = xHatBorio;
        GPSSignals(3,:) = mixtureSignal(:,JNRIndex,nbitsIndex).';
        
        for i = 1:3
            sspace = 0;                 % Search space were the results will be stored
            for ii = 1:K
                y = GPSSignals(i, (ii - 1) * Nc + (1:Nc) );   % use just 1 period of code at the time
                % Compute the search space for a single coherent integration epoch
                Tsspace = DftParallelCodePhaseAcquisition( y, locC, Nc, Nd, DopStep, fs, fi );
                sspace = sspace + Tsspace;  % Non-coherently accumulate the results
            end           
            
            figure
            surf( codeDl*1e6, ((1:Nd) - ceil(Nd/2))*DopStep, sspace ./ max(sspace(:)), 'EdgeColor', 'none');
            axis tight
            xlabel('Delay [$\mu$s]')
            ylabel('Doppler Frequency (Hz)')
            zlabel('Normalized Correlation')
%             if i == 1
%                 title(['NMF ' num2str(JNRVector(JNRIndex)) ' dB']);
%             elseif i == 2
%                 title(['Pulse Blanking ' num2str(JNRVector(JNRIndex)) ' dB']);
%             else
%                 title(['No Mitigation Technique ' num2str(JNRVector(JNRIndex)) ' dB']);
%             end
            
            figProp = struct('size' , fontsize , 'font' , fontname , 'lineWidth' , linewidth, 'figDim', [1 1 1280 720]);
%             formatFig(gcf, [dataPath  'acq_corr_' num2str(i) '_' num2str(JNRVector(JNRIndex))], 'en', figProp);
            
            
        end
    end
    close all;
end

sspace = 0;
for ii = 1:K
    y = locC((ii - 1) * Nc + (1:Nc) );   % use just 1 period of code at the time
    % Compute the search space for a single coherent integration epoch
    Tsspace = DftParallelCodePhaseAcquisition( y, locC, Nc, Nd, DopStep, fs, fi );
    sspace = sspace + Tsspace;  % Non-coherently accumulate the results
end

figure
surf( codeDl*1e6, ((1:Nd) - ceil(Nd/2))*DopStep, sspace./max(sspace(:)), 'EdgeColor', 'none');
axis tight
xlabel('Delay [$\mu$s]')
ylabel('Doppler Frequency (Hz)')
zlabel('Normalized Correlation')
title(['GPS Signal']);

figure;
plot(JNRVector, zeroedSamples)
xlabel('JNR [dB]')
ylabel('Discarded samples');
figProp = struct('size' , fontsize , 'font' , fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
% formatFig(gcf, [dataPath  'pb_discarded_samples'], 'en', figProp);


rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'Misc'])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);