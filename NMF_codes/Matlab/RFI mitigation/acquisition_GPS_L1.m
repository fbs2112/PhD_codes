clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');

addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'Misc'])
addpath(['..' filesep 'Sigtools' filesep])

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
dataPath = ['..' filesep 'figs' filesep '2020-10-05' filesep];
figProp = struct('size' , fontsize , 'font' , fontname , 'lineWidth' , linewidth, 'figDim', [1 1 1280 720]);

load(['..' filesep  'signalsGeneration' filesep 'sim_params' filesep 'sim_params_3.mat']);

matObj = matfile(['.' filesep 'data' filesep 'nmf_testing_3.mat']);

fi = 0;
fs = paramsSignal.Freqsamp;
Nd = 81;
DopStep = 125; %(Hz)

trueDelay = 10e-6;
trueDoppler = 1.2e3;
paramsSignal.numberOfGPSSignals = 1;

[~,locC] = GPSGen(paramsSignal);   % Resample the code at data sampling frequency

Tc = paramsSignal.Intetime;     % Coherent integration time in ms
Nc = floor(Tc * fs);

K = 1;  % Number of non-coherent integrations
codeDl = (0:(Nc - 1)) / fs; % Code delays

if bitand( Nd, 1 ) == 1    % It is an odd number
    Freq = (-((Nd - 1) / 2):((Nd - 1) / 2) ) * DopStep;
else
    Freq = (-(Nd/2):( (Nd-2) / 2 ) ) * DopStep;
end

loopIndex = 1;
bandwidthVector = (2:6:14);
JNRVector = matObj.JNRVector;
for bandwidthIndex = 1:length(bandwidthVector)
    for JNRIndex = 1:length(JNRVector)
        
        GPSSignals(1,:) = matObj.xHat(:,2,JNRIndex,1,bandwidthIndex,loopIndex).';
%         GPSSignals(2,:) = matObj.xHat2(:,1,JNRIndex,1,bandwidthIndex,loopIndex).';
        
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
            title(['Bandwidth: ' num2str(bandwidthVector(bandwidthIndex)) ' MHz' ' JNR: ' num2str(JNRVector(JNRIndex)) ' dB'])
            
            [maxVal, DopInd] = max(max(sspace.'));
            DopFreq(JNRIndex,i) = Freq(DopInd);
            [maxVal, codInd] = max(max(sspace));
            delay(JNRIndex,i) = codeDl(codInd)*1e6;
        end
        close all;
        
    end
    fprintf('-------------------------\n');
    
end
% sspace = 0;
% locCDelayed = [locC(end - round(trueDelay*fs)+1:end) locC(1:end - round(trueDelay*fs))];
% Timeofthisloop = 0:Nc-1;
% Carrphase = mod(2*pi*(2e3)*Timeofthisloop/fs,2*pi);
% Carrier = exp(1i*Carrphase);
% locCDelayed = locCDelayed.*Carrier;
%
% for ii = 1:K
%     y = locCDelayed((ii - 1) * Nc + (1:Nc) );   % use just 1 period of code at the time
%     % Compute the search space for a single coherent integration epoch
%     Tsspace = DftParallelCodePhaseAcquisition( y, locC, Nc, Nd, DopStep, fs, fi );
%     sspace = sspace + Tsspace;  % Non-coherently accumulate the results
% end
%
% figProp = struct('size' , fontsize , 'font' , fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
%
% figure
% surf( codeDl*1e6, ((1:Nd) - ceil(Nd/2))*DopStep/1e3, sspace./max(sspace(:)), 'EdgeColor', 'none');
% axis tight
% xlabel('Delay [$\mu$s]')
% ylabel('Doppler Frequency (kHz)')
% zlabel('Normalised Correlation')
% formatFig(gcf, [dataPath 'caf'], 'en', figProp);


rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'Misc'])
rmpath(['..' filesep 'signalsGeneration' filesep]);