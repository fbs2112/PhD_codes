clear;
clc;
close all;

addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'Sigtools' filesep])

load(['..' filesep  'signalsGeneration' filesep 'sim_params' filesep 'sim_params_3.mat']);

matObj1 = matfile(['.' filesep 'data' filesep 'nmf_testing_15.mat']);                      %stft
matObj2 = matfile(['.' filesep 'data' filesep 'nmf_testing_16.mat']);                      %fsst
matObj3 = matfile(['.' filesep 'data' filesep 'nmf_testing_17.mat']);                      %stft semi
matObj4 = matfile(['.' filesep 'data' filesep 'nmf_testing_18.mat']);                      %fsst semi
matObj5 = matfile(['..' filesep 'pulsedInterf' filesep 'data' filesep 'resultsPai11.mat']);%kalman
matObj6 = matfile(['..' filesep 'pulsedInterf' filesep 'data' filesep 'resultsPai10.mat']);%notch
matObj7 = matfile(['.' filesep 'data' filesep 'wav_testing_2.mat']);                       %wavelet
matObj8 = matfile(['.' filesep 'data' filesep 'nmf_testing_19.mat']);
matObj9 = matfile(['.' filesep 'data' filesep 'nmf_testing_20.mat']);

fi = 0;
fs = paramsSignal.Freqsamp;
Nd = 81;
DopStep = 125; %(Hz)

trueDelay = 10e-6;
trueDoppler = 1e3;
[~,locC] = GPSGen(paramsSignal);
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

bandwidthVector = [2 8 14];
monteCarloLoops = 100;
JNRVector = matObj3.JNRVector;

corrOut = zeros(11, length(JNRVector), length(bandwidthVector), monteCarloLoops);
generalisedSNR = zeros(11, length(JNRVector), length(bandwidthVector), monteCarloLoops);
peakRatio = zeros(11, length(JNRVector), length(bandwidthVector), monteCarloLoops);
DopFreq = zeros(11, length(JNRVector), length(bandwidthVector), monteCarloLoops);
codeDelay = zeros(11, length(JNRVector), length(bandwidthVector), monteCarloLoops);

for loopIndex = 1:monteCarloLoops
    loopIndex
    xHatAux = squeeze(matObj1.xHat(:,2,:,1,:,loopIndex));
    xHatAux2 = squeeze(matObj2.xHat(:,2,:,1,:,loopIndex));
    xHatAux3 = squeeze(matObj3.xHat(:,2,:,1,:,loopIndex));
    xHatAux4 = squeeze(matObj4.xHat(:,2,:,1,:,loopIndex));
    xHatPaiPerfAuxKalman = squeeze(matObj5.xHatPaiPerf(:,:,1,1,:,loopIndex));
    xHatPaiIFestAuxKalman = squeeze(matObj5.xHatPaiIFest(:,:,1,1,:,loopIndex));
    xHatPaiPerfAuxNotch = squeeze(matObj6.xHatPaiPerf(:,:,1,1,:,loopIndex));
    xHatPaiIFestAuxNotch = squeeze(matObj6.xHatPaiIFest(:,:,1,1,:,loopIndex));
    xHatAux5 = squeeze(matObj7.xHat(:,:,1,:,loopIndex));
    xHatAux6 = squeeze(matObj8.xHat(:,2,:,1,:,loopIndex));
    xHatAux7 = squeeze(matObj9.xHat(:,2,:,1,:,loopIndex));
    
    for bandwidthIndex = 1:length(bandwidthVector)
        bandwidthIndex
        for JNRIndex = 1:length(JNRVector)
            
            GPSSignals(1,:) = xHatAux(:,JNRIndex,bandwidthIndex);
            GPSSignals(2,:) = xHatAux2(:,JNRIndex,bandwidthIndex);
            GPSSignals(3,:) = xHatAux3(:,JNRIndex,bandwidthIndex);
            GPSSignals(4,:) = xHatAux4(:,JNRIndex,bandwidthIndex);
            GPSSignals(5,:) = xHatPaiPerfAuxKalman(:,JNRIndex,bandwidthIndex);
            GPSSignals(6,:) = xHatPaiIFestAuxKalman(:,JNRIndex,bandwidthIndex);
            GPSSignals(7,:) = xHatPaiPerfAuxNotch(:,JNRIndex,bandwidthIndex);
            GPSSignals(8,:) = xHatPaiIFestAuxNotch(:,JNRIndex,bandwidthIndex);
            GPSSignals(9,:) = xHatAux5(:,JNRIndex,bandwidthIndex);
            GPSSignals(10,:) = xHatAux6(:,JNRIndex,bandwidthIndex);
            GPSSignals(11,:) = xHatAux7(:,JNRIndex,bandwidthIndex);
            for i = 1:11
                sspace = 0;                 % Search space were the results will be stored
                for ii = 1:K
                    y = GPSSignals(i, (ii - 1) * Nc + (1:Nc) );   % use just 1 period of code at the time
                    % Compute the search space for a single coherent integration epoch
                    [Tsspace, searchSpace] = DftParallelCodePhaseAcquisition( y, locC, Nc, Nd, DopStep, fs, fi );
                    sspace = sspace + Tsspace;  % Non-coherently accumulate the results
                end
                
                [maxVal, DopInd] = max(max(sspace.'));
                [~, codInd] = max(max(sspace));
                
                corrOut(i,JNRIndex,bandwidthIndex,loopIndex) = searchSpace(DopInd,codInd);
                
                correctPeak = sspace(idxFreq,idxDelay);
                sspaceAux = sspace(:);
                sspaceAux = sspaceAux(sspaceAux~=correctPeak);
                
                peakAux = (mean(correctPeak - sspaceAux)^2);
                varSspace = var(sspace(:));
                generalisedSNR(i,JNRIndex,bandwidthIndex,loopIndex) = peakAux/varSspace;
                
                 %--- Find 1 chip wide C/A code phase exclude range around the peak ----
                samplesPerCodeChip = round(fs / paramsSignal.Freqcode);
                excludeRangeIndex1 = codInd - samplesPerCodeChip;
                excludeRangeIndex2 = codInd + samplesPerCodeChip;
                
                %--- Correct C/A code phase exclude range if the range includes array
                %boundaries
                if excludeRangeIndex1 < 2
                    codePhaseRange = excludeRangeIndex2 : ...
                        (Nc + excludeRangeIndex1);
                    
                elseif excludeRangeIndex2 >= Nc
                    codePhaseRange = (excludeRangeIndex2 - Nc) : ...
                        excludeRangeIndex1;
                else
                    codePhaseRange = [1:excludeRangeIndex1, ...
                        excludeRangeIndex2 : Nc];
                end
                
                %--- Find the second highest correlation peak in the same freq. bin ---
                secondPeakSize = max(sspace(DopInd, codePhaseRange));
                
                peakRatio(i,JNRIndex,bandwidthIndex,loopIndex) = maxVal / secondPeakSize;
                DopFreq(i,JNRIndex,bandwidthIndex,loopIndex) = Freq(DopInd);
                codeDelay(i,JNRIndex,bandwidthIndex,loopIndex) = codeDl(codInd)*1e6;
            end
        end
    end
end

save(['.' filesep 'data' filesep 'fig_merit_6.mat'], 'corrOut', 'generalisedSNR', 'peakRatio', 'DopFreq', 'codeDelay', '-v7.3');

rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep]);