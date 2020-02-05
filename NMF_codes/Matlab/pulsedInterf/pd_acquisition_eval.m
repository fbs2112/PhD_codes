addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
addpath(['..' filesep 'Misc'])
addpath(['..' filesep 'Sigtools' filesep])

%%
%Probability of detection evaluation
clear;
clc;
close all;
    
load sim_params_2.mat;
load(['.' filesep 'data' filesep 'nmf_testing_16.mat']);

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
% JNRVector = [0];

Tc = paramsSignal.Intetime;     % Coherent integration time in ms
Nc = floor(Tc * fs);

K = 1;  % Number of non-coherent integrations
codeDl = (0:(Nc - 1)) / fs; % Code delays
Carrphase = mod(2*pi*(trueDoppler)*codeDl,2*pi);
Carrier = exp(1i*Carrphase).';
monteCarloLoops = size(mixtureSignal, 3);

thresPeakRatio = linspace(1e5, 9e5, 50);
thresPeakRatioBorre = linspace(1, 10, 50);
thresPAPR = linspace(5, 80, 50);

% for nbitsIndex = 1:length(nbits)
for loopIndex = 1:monteCarloLoops
    loopIndex
%     locC = quantise_gps(locC, nbits(nbitsIndex));
    for JNRIndex = 1:length(JNR)
        
        xHatBorio = mixtureSignal(:,JNRIndex,loopIndex);
        zeroedSamples(JNRIndex,loopIndex) = length(find(abs(xHatBorio) > thresholdBorio));
        xHatBorio(abs(xHatBorio) > thresholdBorio) = 0;
        
        GPSSignals(1,:) = xHat(:,2,JNRIndex,loopIndex).';
        GPSSignals(2,:) = xHatBorio;
        GPSSignals(3,:) = mixtureSignal(:,JNRIndex,loopIndex).';
        
        for i = 1:3
            sspace = 0;                 % Search space were the results will be stored
            for ii = 1:K
                y = GPSSignals(i, (ii - 1) * Nc + (1:Nc) );   % use just 1 period of code at the time
                % Compute the search space for a single coherent integration epoch
                Tsspace = DftParallelCodePhaseAcquisition( y, locC, Nc, Nd, DopStep, fs, fi );
                sspace = sspace + Tsspace;  % Non-coherently accumulate the results
            end           
            signalPower = GPSSignals(i,:)*GPSSignals(i,:)'/length(GPSSignals(i,:));
            [maxVal(i,JNRIndex,loopIndex), DopInd] = max(max(sspace.'));
            peakRatio(i,JNRIndex,loopIndex) = maxVal(i,JNRIndex,loopIndex)/signalPower;
            [peakVector, idx] = sort(sspace(DopInd,:), 'descend');
            idxAux = idx(2:end);
            minPeakDistance = round(fs*1e-6);
            secondPeakIndexAux = find(abs(idxAux - idx(1)) > minPeakDistance, 1, 'first');
            peakRatioBorre(i,JNRIndex,loopIndex) = maxVal(i,JNRIndex,loopIndex)/sspace(DopInd,idxAux(secondPeakIndexAux));
            
            papr(i,JNRIndex,loopIndex) = maxVal(i,JNRIndex,loopIndex) / mean(sspace(:));
            
            
            GPSSignalNoDoppler = (GPSSignals(i,:).').*conj(Carrier);
            locCDelayed = [locC(end - round(trueDelay*fs)+1:end) locC(1:end - round(trueDelay*fs))]; % Introducing artificial code delay
            navSignalHat(i,JNRIndex,loopIndex) = GPSSignalNoDoppler.'*(locCDelayed).';
            
        end
    end
end

averagePeakValue = mean(maxVal, 3);
averagePeakRatio = mean(peakRatio, 3);
averagePeakRatioBorre = mean(peakRatioBorre, 3);
averagePapr = mean(papr, 3);

stdPeakValue = std(maxVal, 0, 3);
stdPeakRatio = std(peakRatio, 0, 3);
stdPeakRatioBorre = std(peakRatioBorre, 0, 3);
stdPapr = std(papr, 0, 3);


for i = 1:length(thresPeakRatio)
    tpPR(:,:,:,i) = peakRatio > thresPeakRatio(i);
    fnPR(:,:,:,i) = not(tpPR(:,:,:,i));
    
    tpPRB(:,:,:,i) = peakRatioBorre > thresPeakRatioBorre(i);
    fnPRB(:,:,:,i) = not(tpPRB(:,:,:,i));
    
    tpPAPR(:,:,:,i) = papr > thresPAPR(i);
    fnPAPR(:,:,:,i) = not(tpPAPR(:,:,:,i));
end

for JNRIndex = 1:length(JNR)
    for i = 1:3
        
%         SNRHat(i,JNRIndex) = SNResti_pai(squeeze(navSignalHat(i,JNRIndex,:)));
%         SNRHat_SNV(i,JNRIndex) = SNResti_SNV(squeeze(navSignalHat(i,JNRIndex,:)));
%         figure;
%         histogram(peakRatio(i,JNRIndex,:), 10);
%         figure;
%         histogram(peakRatioBorre(i,JNRIndex,:), 10);
%         figure;
%         histogram(papr(i,JNRIndex,:), 10);
    end
end

% save(['.' filesep 'data' filesep 'detection01.mat'], 'tpPR', 'fnPR', 'tpPRB', 'fnPRB', 'tpPAPR', 'fnPAPR', 'SNRHat', 'SNRHat_SNV');
save(['.' filesep 'data' filesep 'detection01.mat'], 'tpPR', 'fnPR', 'tpPRB', 'fnPRB', 'tpPAPR', 'fnPAPR');



%%
%Probability of false alarm evaluation
clear;
clc;
close all;

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

JNR = [0];

Tc = paramsSignal.Intetime;     % Coherent integration time in ms
Nc = floor(Tc * fs);

K = 1;  % Number of non-coherent integrations
codeDl = (0:(Nc - 1)) / fs; % Code delays

monteCarloLoops = size(mixtureSignal, 2);

thresPeakRatio = linspace(1e5, 9e5, 50);
thresPeakRatioBorre = linspace(1, 10, 50);
thresPAPR = linspace(5, 80, 50);

% for nbitsIndex = 1:length(nbits)
for loopIndex = 1:monteCarloLoops
    loopIndex
%     locC = quantise_gps(locC, nbits(nbitsIndex));
    for JNRIndex = 1:length(JNR)
        
        xHatBorio = mixtureSignal(:,loopIndex);
        zeroedSamples(JNRIndex,loopIndex) = length(find(abs(xHatBorio) > thresholdBorio));
        xHatBorio(abs(xHatBorio) > thresholdBorio) = 0;
        
        GPSSignals(1,:) = xHat(:,2,JNRIndex,loopIndex).';
        GPSSignals(2,:) = xHatBorio;
        GPSSignals(3,:) = mixtureSignal(:,loopIndex).';
        
        for i = 1:3
            sspace = 0;                 % Search space were the results will be stored
            for ii = 1:K
                y = GPSSignals(i, (ii - 1) * Nc + (1:Nc) );   % use just 1 period of code at the time
                % Compute the search space for a single coherent integration epoch
                Tsspace = DftParallelCodePhaseAcquisition( y, locC, Nc, Nd, DopStep, fs, fi );
                sspace = sspace + Tsspace;  % Non-coherently accumulate the results
            end           
            signalPower = GPSSignals(i,:)*GPSSignals(i,:)'/length(GPSSignals(i,:));
            [maxVal(i,JNRIndex,loopIndex), DopInd] = max(max(sspace.'));
            peakRatio(i,JNRIndex,loopIndex) = maxVal(i,JNRIndex,loopIndex)/signalPower;
            [peakVector, idx] = sort(sspace(DopInd,:), 'descend');
            idxAux = idx(2:end);
            minPeakDistance = round(fs*1e-6);
            secondPeakIndexAux = find(abs(idxAux - idx(1)) > minPeakDistance, 1, 'first');
            peakRatioBorre(i,JNRIndex,loopIndex) = maxVal(i,JNRIndex,loopIndex)/sspace(DopInd,idxAux(secondPeakIndexAux));
            
            papr(i,JNRIndex,loopIndex) = maxVal(i,JNRIndex,loopIndex) / mean(sspace(:));
        end
    end
end
averagePeakValue = mean(maxVal, 3);
averagePeakRatio = mean(peakRatio, 3);
averagePeakRatioBorre = mean(peakRatioBorre, 3);
averagePapr = mean(papr, 3);

stdPeakValue = std(maxVal, 0, 3);
stdPeakRatio = std(peakRatio, 0, 3);
stdPeakRatioBorre = std(peakRatioBorre, 0, 3);
stdPapr = std(papr, 0, 3);


for i = 1:length(thresPeakRatio)
    fpPR(:,:,:,i) = peakRatio > thresPeakRatio(i);
    tnPR(:,:,:,i) = not(fpPR(:,:,:,i));
    
    fpPRB(:,:,:,i) = peakRatioBorre > thresPeakRatioBorre(i);
    tnPRB(:,:,:,i) = not(fpPRB(:,:,:,i));
    
    fpPAPR(:,:,:,i) = papr > thresPAPR(i);
    tnPAPR(:,:,:,i) = not(fpPAPR(:,:,:,i));
end


% for JNRIndex = 1:length(JNR)
%     for i = 1:3
%         figure;
%         histogram(peakRatio(i,JNRIndex,:), 10);
%         figure;
%         histogram(peakRatioBorre(i,JNRIndex,:), 10);
%         figure;
%         histogram(papr(i,JNRIndex,:), 10);
%     end
% end


load(['.' filesep 'data' filesep 'detection01.mat']);
% save(['.' filesep 'data' filesep 'detection01.mat'], 'tpPR', 'fnPR', 'tpPRB', 'fnPRB', 'tpPAPR', 'fnPAPR', ...
%     'fpPR', 'tnPR', 'fpPRB', 'tnPRB', 'fpPAPR', 'tnPAPR', 'SNRHat', 'SNRHat_SNV');
save(['.' filesep 'data' filesep 'detection01.mat'], 'tpPR', 'fnPR', 'tpPRB', 'fnPRB', 'tpPAPR', 'fnPAPR', ...
    'fpPR', 'tnPR', 'fpPRB', 'tnPRB', 'fpPAPR', 'tnPAPR');


%%

rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'Misc'])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);