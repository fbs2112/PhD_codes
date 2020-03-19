addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
addpath(['..' filesep 'Misc'])
addpath(['..' filesep 'Sigtools' filesep])

%Probability of detection evaluation for quantised signals
clear;
clc;
close all;

load sim_params_2.mat;
load(['.' filesep 'data' filesep 'nmf_testing_19.mat']);

% varN = 2;
% Pfa = 0.05;
% thresholdBorio = sqrt(-2*varN*log(Pfa));
thresholdBorio = [2 5];
fi = 0;
fs = paramsSignal.Freqsamp;
Nd = 81;
DopStep = 125; %(Hz)

trueDelay = 10e-6;
trueDoppler = 1.5e3;
paramsSignal.numberOfGPSSignals = 1;

[~, locC] = GPSGen(paramsSignal);   % Resample the code at data sampling frequency

Tc = paramsSignal.Intetime;     % Coherent integration time in ms
Nc = floor(Tc * fs);

K = 1;  % Number of non-coherent integrations
codeDl = (0:(Nc - 1)) / fs; % Code delays
Carrphase = mod(2*pi*(trueDoppler)*codeDl,2*pi);
Carrier = exp(1i*Carrphase).';
monteCarloLoops = size(mixtureSignal, 4);

thresPeakRatio = linspace(1e5, 9e5, 50);
thresPeakRatioBorre = linspace(1, 10, 50);
thresPAPR = linspace(5, 80, 50);
JNR = JNRVector;
% for nbitsIndex = 1:length(nbits)
for loopIndex = 1:monteCarloLoops
    loopIndex
    for nbitsIndex = 1:length(nbits)
        %     locC = quantise_gps(locC, nbits(nbitsIndex));
        for JNRIndex = 1:length(JNR)
            
            xHatBorio = mixtureSignal(:,JNRIndex,nbitsIndex,loopIndex);
            zeroedSamples(JNRIndex,nbitsIndex,loopIndex) = length(find(abs(xHatBorio) > thresholdBorio(nbitsIndex)));
            xHatBorio(abs(xHatBorio) > thresholdBorio(nbitsIndex)) = 0;
            
            GPSSignals(1,:) = xHat(:,2,JNRIndex,nbitsIndex,loopIndex).';
            GPSSignals(2,:) = xHatBorio;
            GPSSignals(3,:) = mixtureSignal(:,JNRIndex,nbitsIndex,loopIndex).';
            
            for i = 1:3
                sspace = 0;                 % Search space were the results will be stored
                for ii = 1:K
                    y = GPSSignals(i, (ii - 1) * Nc + (1:Nc) );   % use just 1 period of code at the time
                    % Compute the search space for a single coherent integration epoch
                    [Tsspace, TsspaceComplex(i,JNRIndex,nbitsIndex,loopIndex)] = DftParallelCodePhaseAcquisition( y, locC, Nc, Nd, DopStep, fs, fi );
                    sspace = sspace + Tsspace;  % Non-coherently accumulate the results
                end
                signalPower = GPSSignals(i,:)*GPSSignals(i,:)'/length(GPSSignals(i,:));
                [maxVal(i,JNRIndex,nbitsIndex,loopIndex), DopInd] = max(max(sspace.'));
                peakRatio(i,JNRIndex,nbitsIndex,loopIndex) = maxVal(i,JNRIndex,nbitsIndex,loopIndex)/signalPower;
                [peakVector, idx] = sort(sspace(DopInd,:), 'descend');
                idxAux = idx(2:end);
                minPeakDistance = round(fs*1e-6);
                secondPeakIndexAux = find(abs(idxAux - idx(1)) > minPeakDistance, 1, 'first');
                peakRatioBorre(i,JNRIndex,nbitsIndex,loopIndex) = maxVal(i,JNRIndex,nbitsIndex,loopIndex)/sspace(DopInd,idxAux(secondPeakIndexAux));
                
                papr(i,JNRIndex,nbitsIndex,loopIndex) = maxVal(i,JNRIndex,nbitsIndex,loopIndex) / mean(sspace(:));
                
                GPSSignalNoDoppler = (GPSSignals(i,:).').*conj(Carrier);
                locCDelayed = [locC(end - round(trueDelay*fs)+1:end) locC(1:end - round(trueDelay*fs))]; % Introducing artificial code delay
                navSignalHat(i,JNRIndex,nbitsIndex,loopIndex) = (GPSSignalNoDoppler.'*(locCDelayed).');
                
            end
        end
    end
end

averagePeakValue = mean(maxVal, 4);
averagePeakRatio = mean(peakRatio, 4);
averagePeakRatioBorre = mean(peakRatioBorre, 4);
averagePapr = mean(papr, 4);

stdPeakValue = std(maxVal, 0, 4);
stdPeakRatio = std(peakRatio, 0, 4);
stdPeakRatioBorre = std(peakRatioBorre, 0, 4);
stdPapr = std(papr, 0, 4);


% for i = 1:length(thresPeakRatio)
%     tpPR(:,:,:,i) = peakRatio > thresPeakRatio(i);
%     fnPR(:,:,:,i) = not(tpPR(:,:,:,i));
%     
%     tpPRB(:,:,:,i) = peakRatioBorre > thresPeakRatioBorre(i);
%     fnPRB(:,:,:,i) = not(tpPRB(:,:,:,i));
%     
%     tpPAPR(:,:,:,i) = papr > thresPAPR(i);
%     fnPAPR(:,:,:,i) = not(tpPAPR(:,:,:,i));
% end

for nbitsIndex = 1:length(nbits)
    for JNRIndex = 1:length(JNR)
        for i = 1:3

            SNRHat(i,JNRIndex,nbitsIndex) = SNResti_pai(squeeze(navSignalHat(i,JNRIndex,nbitsIndex,:)));
            SNRHat_SNV(i,JNRIndex,nbitsIndex) = SNResti_SNV(squeeze(navSignalHat(i,JNRIndex,nbitsIndex,:)));
%             SNRHat_MM(i,JNRIndex,nbitsIndex) = SNResti_MM(squeeze(navSignalHat(i,JNRIndex,nbitsIndex,:)));
            SNRHat_NWPR(i,JNRIndex,nbitsIndex) = SNResti_NWPR(squeeze(navSignalHat(i,JNRIndex,nbitsIndex,:)), 10);
            %         figure;
            %         histogram(peakRatio(i,JNRIndex,:), 10);
            %         figure;
            %         histogram(peakRatioBorre(i,JNRIndex,:), 10);
            %         figure;
            %         histogram(papr(i,JNRIndex,:), 10);
        end
    end
end

% save(['.' filesep 'data' filesep 'detection01.mat'], 'tpPR', 'fnPR', 'tpPRB', 'fnPRB', 'tpPAPR', 'fnPAPR', 'SNRHat', 'SNRHat_SNV');
% save(['.' filesep 'data' filesep 'detection01.mat'], 'tpPR', 'fnPR', 'tpPRB', 'fnPRB', 'tpPAPR', 'fnPAPR');







rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'Misc'])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);