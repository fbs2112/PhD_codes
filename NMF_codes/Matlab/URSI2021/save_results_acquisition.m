clear;
clc;
close all;

addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
addpath(['..' filesep 'Misc'])
addpath(['..' filesep 'Sigtools' filesep])

load(['.' filesep 'GPS signals' filesep 'GPS_L5_train']);
load sim_params_L5_1.mat;
load(['.' filesep 'data' filesep 'nmf_testing_11_merged.mat']);
xOut1024 = xOut;
mixtureSignalMerged1024 = mixtureSignalMerged;
load(['.' filesep 'data' filesep 'nmf_testing_10_merged.mat']);

fi = 0;
fs = paramsSignal.Freqsamp;
Nd = 81;
DopStep = 125; %(Hz)

varNoise = 1;
pfa = 0.01;
PB_threshold = sqrt(-2*varNoise*log(pfa));

trueDelay = 10e-6;
trueDoppler = 1.5e3;

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

GPSL5 = GPSSignals(1:Nc);

monteCarloLoops = 200;
JNRVector = JNRVector + 20;

totalSamples = paramsSignal.Intetime * fs;
Timeofthisloop = 0:totalSamples-1;

numberOfZeroedSamples = zeros(length(JNRVector), monteCarloLoops);
numberOfZeroedSamplesFreq = zeros(length(JNRVector), monteCarloLoops);
maxVal = zeros(5, length(JNRVector), monteCarloLoops);
DopInd = zeros(5, length(JNRVector), monteCarloLoops);
codInd = zeros(5, length(JNRVector), monteCarloLoops);
corrOut = zeros(5, length(JNRVector), monteCarloLoops);
generalisedSNR = zeros(5, length(JNRVector), monteCarloLoops);
generalisedSNRMaxPeak  = zeros(5, length(JNRVector), monteCarloLoops);
PPR = zeros(5, length(JNRVector), monteCarloLoops);
peakRatioBorre = zeros(5, length(JNRVector), monteCarloLoops);

for loopIndex = 1:monteCarloLoops
    loopIndex
    for JNRIndex = 1:length(JNRVector)
        
        xAux(:,1) = xOut(:,1,JNRIndex,1,loopIndex);
        
        aux = mixtureSignalMerged(:,JNRIndex,loopIndex);
        numberOfZeroedSamples(JNRIndex,loopIndex) = length(find(abs(aux) > PB_threshold));
        aux(abs(aux) > PB_threshold) = 0;
        xAux(:,3) = aux;
        aux = fft(mixtureSignalMerged(:,JNRIndex,loopIndex));
        numberOfZeroedSamplesFreq(JNRIndex,loopIndex) = length(find((abs(aux)/sqrt(length(aux))) > PB_threshold));
        aux((abs(aux)/sqrt(length(aux))) > PB_threshold) = 0;
        xAux(:,4) = ifft(aux);
        xAux(:,5) = mixtureSignalMerged(:,JNRIndex,loopIndex);
        xAux(:,2) = xOut1024(:,1,JNRIndex,1,loopIndex);
        for i = 1:5
            sspace = 0;                 % Search space were the results will be stored
            for ii = 1:K
                y = xAux((ii - 1) * Nc + (1:Nc),i);   % use just 1 period of code at the time
                % Compute the search space for a single coherent integration epoch
                [Tsspace, searchSpace] = DftParallelCodePhaseAcquisition( y.', GPSL5.', Nc, Nd, DopStep, fs, fi );
                sspace = sspace + Tsspace;  % Non-coherently accumulate the results
            end
            
            [maxVal(i,JNRIndex,loopIndex), DopInd(i,JNRIndex,loopIndex)] = max(max(sspace.'));
            [~, codInd(i,JNRIndex,loopIndex)] = max(max(sspace));
            
            corrOut(i,JNRIndex,loopIndex) = searchSpace(DopInd(i,JNRIndex,loopIndex),codInd(i,JNRIndex,loopIndex));
                        
            correctPeak = sspace(idxFreq,idxDelay);
            
            maxPeak = sspace(DopInd(i,JNRIndex,loopIndex),codInd(i,JNRIndex,loopIndex));
            
            sspaceAux = sspace(:);
            sspaceAux = sspaceAux(sspaceAux~=correctPeak);
            
            sspaceAux2 = sspace(:);
            sspaceAux2 = sspaceAux2(sspaceAux2~=maxPeak);
                                   
            peakAux = (mean(correctPeak - sspaceAux)^2);
            varSspace = var(sspace(:));
            
            peakAux2 = (mean(maxPeak - sspaceAux2)^2);
            
            generalisedSNR(i,JNRIndex,loopIndex) = peakAux/varSspace;
            generalisedSNRMaxPeak(i,JNRIndex,loopIndex) = peakAux2 / varSspace;
            
            [pks,locs] = findpeaks(sspace(:), 'SortStr', 'descend');
            [r,c] = ind2sub(size(sspace), locs);
             
            PPR(i,JNRIndex,loopIndex) = pks(1) / pks(2); 
            
            peakMtx(i,JNRIndex,loopIndex,:) = [pks(1) pks(2)];
            [peakVector, idx] = sort(sspace(DopInd(i,JNRIndex,loopIndex),:), 'descend');
            idxAux = idx(2:end);
            minPeakDistance = round(fs*1e-6);
            secondPeakIndexAux = find(abs(idxAux - idx(1)) > minPeakDistance, 1, 'first');
            peakRatioBorre(i,JNRIndex,loopIndex) = maxVal(i,JNRIndex,loopIndex)/sspace(DopInd(i,JNRIndex,loopIndex),idxAux(secondPeakIndexAux));
            
        end
    end
end

M = 20;
for JNRIndex = 1:length(JNRVector)
    for i = 1:5
        corrBuffer = buffer(corrOut(i,JNRIndex,:), M, 0, 'nodelay');
        Pn = sum(real(corrBuffer), 1).^2 + sum(imag(corrBuffer), 1).^2;
        Pw = sum(real(corrBuffer).^2 + imag(corrBuffer).^2, 1);
        s(i,JNRIndex) = mean(Pn ./ Pw);
        stdS(i,JNRIndex) = std(Pn ./ Pw);
        if s(i,JNRIndex) < 1
            s(i,JNRIndex) = 1;
        end
        C_N0(i,JNRIndex) = 10*log10((1/1e-3) * ((s(i,JNRIndex) - 1))/((M - s(i,JNRIndex))));
    end
end

averagePPR = mean(PPR, 3);
averageGenSNRMaxPeak = mean(generalisedSNRMaxPeak, 3);
averageGenSNR = mean(generalisedSNR, 3);
peakToPeakRatio = mean(peakRatioBorre, 3);

save(['.' filesep 'data' filesep 'results_10_11'], 'averagePPR', 'averageGenSNRMaxPeak', 'averageGenSNR', 'peakToPeakRatio', 's', 'C_N0');

rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'Misc'])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);