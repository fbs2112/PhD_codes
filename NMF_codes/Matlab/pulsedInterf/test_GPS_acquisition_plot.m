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
dataPath = ['..' filesep 'figs' filesep '2020-06-05' filesep];

load sim_params_3.mat;

% load(['.' filesep 'data' filesep 'resultsPai02.mat']);
% xHatPaiPerf = xHatPai;
load(['.' filesep 'data' filesep 'nmf_testing_33_1.mat']);
% load(['.' filesep 'data' filesep 'resultsPai03.mat']);
% xHatPaiEst = xHatPai;

% thresholdBorio = [1.5 5];   %for real case
thresholdBorio = [1.5 11];   %for complex case

fi = 0;
fs = paramsSignal.Freqsamp;
Nd = 81;
DopStep = 125; %(Hz)

trueDelay = 10e-6;
trueDoppler = 1.2e3;
paramsSignal.numberOfGPSSignals = 1;

[v, locCAux] = GPSGen(paramsSignal);   % Resample the code at data sampling frequency

Tc = paramsSignal.Intetime;     % Coherent integration time in ms
Nc = floor(Tc * fs);

K = 1;  % Number of non-coherent integrations
codeDl = (0:(Nc - 1)) / fs; % Code delays

if bitand( Nd, 1 ) == 1    % It is an odd number
    Freq = (-((Nd - 1) / 2):((Nd - 1) / 2) ) * DopStep;
else
    Freq = (-(Nd/2):( (Nd-2) / 2 ) ) * DopStep;
end

nbits = [2];
loopIndex = 1;

% xHat = xHatSemi;
bandwidthVector = (2:3:14)*1e6;
for bandwidthIndex = 1:1%length(bandwidthVector)
    for nbitsIndex = 1:length(nbits)
        locC = quantise_gps(locCAux, nbits(nbitsIndex));
        for JNRIndex = 1:length(JNRVector)
            
            %         xHatBorio = mixtureSignal(:,JNRIndex,nbitsIndex,loopIndex);
            %         zeroedSamples(JNRIndex,nbitsIndex) = length(find(abs(xHatBorio) > thresholdBorio(nbitsIndex)));
            %         fprintf('Number of zeroed samples: %i\n', zeroedSamples(JNRIndex,nbitsIndex));
            %         xHatBorio(abs(xHatBorio) > thresholdBorio(nbitsIndex)) = 0;
            GPSSignals(1,:) = xHat(:,1,JNRIndex,nbitsIndex,bandwidthIndex,loopIndex).';
%             GPSSignals(1,:) = abs(randn(1, Nc) + 1j*randn(1, Nc)) .*exp(1j*v.');
%             GPSSignals(2,:) = xHatPaiPerf(:,JNRIndex,nbitsIndex,bandwidthIndex,loopIndex).';
%             GPSSignals(3,:) = xHatPaiEst(:,JNRIndex,nbitsIndex,bandwidthIndex,loopIndex).';
            %         GPSSignals(2,:) = xHatBorio;
            %         GPSSignals(3,:) = mixtureSignal(:,JNRIndex,nbitsIndex,loopIndex ).';
            
            for i = 1:1
                sspace = 0;                 % Search space were the results will be stored
                for ii = 1:K
                    y = GPSSignals(i, (ii - 1) * Nc + (1:Nc) );   % use just 1 period of code at the time
                    % Compute the search space for a single coherent integration epoch
                    Tsspace = DftParallelCodePhaseAcquisition( y, locC, Nc, Nd, DopStep, fs, fi );
                    sspace = sspace + Tsspace;  % Non-coherently accumulate the results
                end
                
                figure
                surf( codeDl*1e6, ((1:Nd) - ceil(Nd/2))*DopStep, sspace./max(sspace(:)), 'EdgeColor', 'none');
                axis tight
                xlabel('Delay [$\mu$s]')
                ylabel('Doppler Frequency (Hz)')
                zlabel('Normalised Correlation')
               
%                 if i == 2
%                     text(350, 0.8, 1, [num2str(zeroedSamples(JNRIndex,nbitsIndex)*100/Nc, '%2.2f') '\% discarded samples'], 'Fontsize', 14);
%                 end
                figProp = struct('size' , fontsize , 'font' , fontname , 'lineWidth' , linewidth, 'figDim', [1 1 1280 720]);
%                  formatFig(gcf, [dataPath  'acq_corr_' num2str(i) '_' num2str(JNRVector(JNRIndex)) '_' num2str(bandwidthIndex)], 'en', figProp);
                
                [maxVal, DopInd] = max(max(sspace.'));
                DopFreq(JNRIndex,i,nbitsIndex) = Freq(DopInd);
                [maxVal, codInd] = max(max(sspace));
                delay(JNRIndex,i,nbitsIndex) = codeDl(codInd)*1e6;
%                 close all;
            end
            
        end
        fprintf('-------------------------\n');
    end
end
sspace = 0;
locCDelayed = [locC(end - round(trueDelay*fs)+1:end) locC(1:end - round(trueDelay*fs))];
for ii = 1:K
    y = locCDelayed((ii - 1) * Nc + (1:Nc) );   % use just 1 period of code at the time
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
plot(JNRVector, zeroedSamples/Nc * 100)
xlabel('JNR [dB]')
ylabel('Discarded samples [\%]');
figProp = struct('size' , fontsize , 'font' , fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
grid on;
xlim([-5 50]);
ylim([5 35]);
legend('2 bits', '4 bits', 'location', 'best');
formatFig(gcf, [dataPath  'pb_discarded_samples'], 'en', figProp);


rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'Misc'])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);