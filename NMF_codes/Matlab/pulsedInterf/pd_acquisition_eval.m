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
dataPath = ['..' filesep 'figs' filesep '2020-06-10' filesep];
figProp = struct('size' , fontsize , 'font' , fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);

load sim_params_3.mat;

load(['.' filesep 'data' filesep 'resultsPai02.mat']);
xHatPaiPerf = xHatPai;
load(['.' filesep 'data' filesep 'nmf_testing_32.mat']);
load(['.' filesep 'data' filesep 'resultsPai03.mat']);
xHatPaiEst = xHatPai;
xHatEst = xHat;
load(['.' filesep 'data' filesep 'nmf_testing_34.mat']);

% xHat = xHatSemi;
% thresholdBorio = [1.5 5];   %for real case
thresholdBorio = [1.5 11];   %for complex case
fi = 0;
fs = paramsSignal.Freqsamp;
Nd = 81;
DopStep = 125; %(Hz)

trueDelay = 10e-6;
trueDoppler = 1e3;

paramsSignal.numberOfGPSSignals = 1;

[~, locC] = GPSGen(paramsSignal);   % Resample the code at data sampling frequency

Tc = paramsSignal.Intetime;     % Coherent integration time in ms
Nc = floor(Tc * fs);

K = 1;  % Number of non-coherent integrations
codeDl = (0:(Nc - 1)) / fs; % Code delays

if bitand( Nd, 1 ) == 1    % It is an odd number
    Freq = (-((Nd - 1) / 2):((Nd - 1) / 2) ) * DopStep;
else
    Freq = (-(Nd/2):( (Nd-2) / 2 ) ) * DopStep;
end

[~,idxFreq] = min(abs(Freq - trueDoppler));
[~,idxDelay] = min(abs(codeDl - trueDelay));

Carrphase = mod(2*pi*(trueDoppler)*codeDl,2*pi);
Carrier = exp(1i*Carrphase).';
monteCarloLoops = 100;

thresPeakRatio = linspace(1e5, 9e5, 50);
thresPeakRatioBorre = linspace(1, 10, 50);
thresPAPR = linspace(5, 80, 50);
JNR = JNRVector;
nbits = 2;

bandwidthVector = (2:3:14)*1e6;

for bandwidthIndex = 1:length(bandwidthVector)
    for loopIndex = 1:monteCarloLoops
        loopIndex
        for nbitsIndex = 1:length(nbits)
            locC = quantise_gps(locC, nbits(nbitsIndex));
            for JNRIndex = 1:length(JNR)
                
                %             xHatBorio = mixtureSignal(:,JNRIndex,nbitsIndex,loopIndex);
                %             zeroedSamples(JNRIndex,nbitsIndex,loopIndex) = length(find(abs(xHatBorio) > thresholdBorio(nbitsIndex)));
                %             xHatBorio(abs(xHatBorio) > thresholdBorio(nbitsIndex)) = 0;
                
                %             GPSSignals(1,:) = xHat(:,JNRIndex,nbitsIndex,loopIndex).';
                GPSSignals(1,:) = xHat(:,2,JNRIndex,nbitsIndex,bandwidthIndex,loopIndex).';
                GPSSignals(2,:) = xHatEst(:,2,JNRIndex,nbitsIndex,bandwidthIndex,loopIndex).';
                GPSSignals(3,:) = xHatPaiPerf(:,JNRIndex,nbitsIndex,bandwidthIndex,loopIndex).';
                GPSSignals(4,:) = xHatPaiEst(:,JNRIndex,nbitsIndex,bandwidthIndex,loopIndex).';
                %             GPSSignals(3,:) = mixtureSignal(:,JNRIndex,nbitsIndex,loopIndex).';
                
                for i = 1:4
                    sspace = 0;                 % Search space were the results will be stored
                    for ii = 1:K
                        y = GPSSignals(i, (ii - 1) * Nc + (1:Nc) );   % use just 1 period of code at the time
                        % Compute the search space for a single coherent integration epoch
                        [Tsspace] = DftParallelCodePhaseAcquisition( y, locC, Nc, Nd, DopStep, fs, fi );
                        sspace = sspace + Tsspace;  % Non-coherently accumulate the results
                    end
                    
                    [maxVal(i,JNRIndex,nbitsIndex,loopIndex), DopInd] = max(max(sspace.'));
                    [~, codInd] = max(max(sspace));
                    
                    correctPeak = sspace(idxFreq,idxDelay);
                    
                    signalPower = GPSSignals(i,:)*GPSSignals(i,:)'/length(GPSSignals(i,:));
                    
                    sspaceAux = sspace(:);
                    sspaceAux = sspaceAux(sspaceAux~=correctPeak);
                    
                    peakAux(i,JNRIndex,nbitsIndex,loopIndex) = (mean(correctPeak - sspaceAux)^2);
                    varSspace(i,JNRIndex,nbitsIndex,loopIndex) = var(sspace(:));
                    
                    generalisedSNR(i,JNRIndex,nbitsIndex,loopIndex) = peakAux(i,JNRIndex,nbitsIndex,loopIndex)/varSspace(i,JNRIndex,nbitsIndex,loopIndex);
                    
                    peakRatio(i,JNRIndex,nbitsIndex,loopIndex) = maxVal(i,JNRIndex,nbitsIndex,loopIndex)/signalPower;
                    
                    [peakVector, idx] = sort(sspace(DopInd,:), 'descend');
                    idxAux = idx(2:end);
                    minPeakDistance = round(fs*1e-6);
                    secondPeakIndexAux = find(abs(idxAux - idx(1)) > minPeakDistance, 1, 'first');
                    peakRatioBorre(i,JNRIndex,nbitsIndex,loopIndex) = maxVal(i,JNRIndex,nbitsIndex,loopIndex)/sspace(DopInd,idxAux(secondPeakIndexAux));
                    
                    papr(i,JNRIndex,nbitsIndex,loopIndex) = maxVal(i,JNRIndex,nbitsIndex,loopIndex) / mean(sspace(:));
                    
%                     GPSSignalNoDoppler = (GPSSignals(i,:).').*conj(Carrier);
%                     locCDelayed = [locC(end - round(trueDelay*fs)+1:end) locC(1:end - round(trueDelay*fs))]; % Introducing artificial code delay
%                     navSignalHat(i,JNRIndex,nbitsIndex,loopIndex) = (GPSSignalNoDoppler.'*(locCDelayed).');
                    
                end
            end
        end
    end
    
    averagePeakValue = mean(maxVal, 4);
    averagePeakRatio = mean(peakRatio, 4);
    averagePeakRatioBorre = mean(peakRatioBorre, 4);
    averagePapr = mean(papr, 4);
    averageGenSNR = mean(generalisedSNR, 4);
    
    stdPeakValue = std(maxVal, 0, 4);
    stdPeakRatio = std(peakRatio, 0, 4);
    stdPeakRatioBorre = std(peakRatioBorre, 0, 4);
    stdPapr = std(papr, 0, 4);
    stdGenSNR = std(generalisedSNR, 0, 4);
    
    figure;
    plot(JNRVector, averagePeakRatioBorre(1,:,1), 'Color', [0 0.4470 0.7410]);
    hold on;
    plot(JNRVector, averagePeakRatioBorre(2,:,1), 'Color', [0.8500 0.3250 0.0980]);
    plot(JNRVector, averagePeakRatioBorre(3,:,1), 'Color', [0.9290 0.6940 0.1250]);
    plot(JNRVector, averagePeakRatioBorre(4,:,1), 'Color', [0.4940, 0.1840, 0.5560]);
    
    grid on;
    xlim([-5 50]);
    ylabel('Borre''s Ratio')
    xlabel('JNR [dB]');
    legend('NMF-based', 'NMF-based est', 'Kalman', 'Kalman-IF est.', 'location', 'best');
%     formatFig(gcf, [dataPath  'av_PRB_' num2str(bandwidthIndex)], 'en', figProp);
    
    figure;
    plot(JNRVector, averagePapr(1,:,1), 'Color', [0 0.4470 0.7410]);
    hold on;
    plot(JNRVector, averagePapr(2,:,1), 'Color', [0.8500 0.3250 0.0980]);
    plot(JNRVector, averagePapr(3,:,1), 'Color', [0.9290 0.6940 0.1250]);
    plot(JNRVector, averagePapr(4,:,1), 'Color', [0.4940, 0.1840, 0.5560]);

    grid on;
    xlim([-5 50]);
    ylabel('PAPR')
    xlabel('JNR [dB]');
    legend('NMF-based', 'NMF-based est', 'Kalman', 'Kalman-IF est.', 'location', 'best');
%     formatFig(gcf, [dataPath  'av_PAPR_' num2str(bandwidthIndex)], 'en', figProp);
    
    figure;
    plot(JNRVector, pow2db(averageGenSNR(1,:,1)), 'Color', [0 0.4470 0.7410]);
    hold on;
    plot(JNRVector, pow2db(averageGenSNR(2,:,1)), 'Color', [0.8500 0.3250 0.0980]);
    plot(JNRVector, pow2db(averageGenSNR(3,:,1)), 'Color', [0.9290 0.6940 0.1250]);
    plot(JNRVector, pow2db(averageGenSNR(4,:,1)), 'Color', [0.4940, 0.1840, 0.5560]);

    grid on;
    xlim([-5 50]);
    ylabel('Generalised SNR [dB]')
    xlabel('JNR [dB]');
    legend('NMF-based', 'NMF-based est', 'Kalman', 'Kalman-IF est.', 'location', 'best');
    formatFig(gcf, [dataPath  'av_SNR_' num2str(bandwidthIndex)], 'en', figProp);
    
    
    % figure;
    % plot(JNRVector, pow2db(averageGenSNR(1,:,1)), 'Color', [0 0.4470 0.7410]);
    % hold on;
    % plot(JNRVector, pow2db(averageGenSNR(1,:,2)), '--', 'Color', [0 0.4470 0.7410]);
    % plot(JNRVector, pow2db(averageGenSNR(2,:,2)), '--', 'Color', [0.8500 0.3250 0.0980]);
    % grid on;
    % xlim([-5 50]);
    % ylim([35 40]);
    % ylabel('Generalised SNR [dB]')
    % xlabel('JNR [dB]');
    % legend('NMF-based, 2 bits', 'NMF-based, 4 bits', 'PB, 4 bits', 'location', 'best');
    % formatFig(gcf, [dataPath  'av_SNR_zoom'], 'en', figProp);
    
    figure;
    histogram(peakRatioBorre(1,4,1,:), 21);
    hold on;
    histogram(peakRatioBorre(2,4,1,:), 21);
    histogram(peakRatioBorre(3,4,1,:), 21);
    histogram(peakRatioBorre(4,4,1,:), 21);
    xlabel('Borre''s Ratio')
    ylabel('Ocurrences');
    legend('NMF-based', 'NMF-based est', 'Kalman', 'Kalman-IF est.');
%     formatFig(gcf, [dataPath  'hist_PRB_' num2str(bandwidthIndex)], 'en', figProp);
    
    figure;
    histogram(papr(1,4,1,:), 21);
    hold on;
    histogram(papr(2,4,1,:), 21);
    histogram(papr(3,4,1,:), 21);
    histogram(papr(4,4,1,:), 21);

    xlabel('PAPR')
    ylabel('Ocurrences');
    legend('NMF-based', 'NMF-based est', 'Kalman', 'Kalman-IF est.');
%     formatFig(gcf, [dataPath  'hist_PAPR_' num2str(bandwidthIndex)], 'en', figProp);
    
    figure;
    histogram(pow2db(generalisedSNR(1,4,1,:)), 21);
    hold on;
    histogram(pow2db(generalisedSNR(2,4,1,:)), 21);
    histogram(pow2db(generalisedSNR(3,4,1,:)), 21);
    histogram(pow2db(generalisedSNR(4,4,1,:)), 21);

    xlabel('Generalised SNR [dB]')
    ylabel('Ocurrences');
    legend('NMF-based', 'NMF-based est', 'Kalman', 'Kalman-IF est.');
    formatFig(gcf, [dataPath  'hist_SNR_' num2str(bandwidthIndex)], 'en', figProp);
    close all;
end

% for nbitsIndex = 1:length(nbits)
%     for JNRIndex = 1:length(JNR)
%         for i = 1:3
%
%             SNRHat(i,JNRIndex,nbitsIndex) = SNResti_pai(squeeze(navSignalHat(i,JNRIndex,nbitsIndex,:)));
%             SNRHat_SNV(i,JNRIndex,nbitsIndex) = SNResti_SNV(squeeze(navSignalHat(i,JNRIndex,nbitsIndex,:)));
%             SNRHat_MM(i,JNRIndex,nbitsIndex) = SNResti_MM(squeeze(navSignalHat(i,JNRIndex,nbitsIndex,:)));
%             SNRHat_NWPR(i,JNRIndex,nbitsIndex) = SNResti_NWPR(squeeze(navSignalHat(i,JNRIndex,nbitsIndex,:)), 10);
%
%         end
%     end
% end

rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'Misc'])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);