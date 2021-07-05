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
dataPath = ['..' filesep 'figs' filesep '2020-06-10' filesep];
figProp = struct('size' , fontsize , 'font' , fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);

load(['..' filesep  'signalsGeneration' filesep 'sim_params' filesep 'sim_params_3.mat']);

matObj1 = matfile(['..' filesep 'pulsedInterf' filesep 'data' filesep 'nmf_testing_35.mat']);
matObj2 = matfile(['..' filesep 'pulsedInterf' filesep 'data' filesep 'resultsPai07.mat']);
matObj3 = matfile(['.' filesep 'data' filesep 'nmf_testing_7.mat']);
matObj4 = matfile(['.' filesep 'data' filesep 'nmf_testing_6.mat']);

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

for loopIndex = 1:monteCarloLoops
    loopIndex
    xHatAux = squeeze(matObj1.xHat(:,2,:,1,1,:,loopIndex));
    xHatPaiPerfAux = squeeze(matObj2.xHatPaiPerf(:,:,1,1,:,loopIndex));
    xHatPaiIFestAux = squeeze(matObj2.xHatPaiIFest(:,:,1,1,:,loopIndex));
    xHatAux2 = squeeze(matObj3.xHat(:,2,:,1,:,loopIndex));
    xHatAux3 = squeeze(matObj4.xHat(:,2,:,1,:,loopIndex));
    for bandwidthIndex = 1:length(bandwidthVector)
        for JNRIndex = 1:length(JNRVector)
            GPSSignals(1,:) = xHatAux(:,JNRIndex,bandwidthIndex);
            GPSSignals(2,:) = xHatPaiPerfAux(:,JNRIndex,bandwidthIndex);
            GPSSignals(3,:) = xHatPaiIFestAux(:,JNRIndex,bandwidthIndex);
            GPSSignals(4,:) = xHatAux2(:,JNRIndex,bandwidthIndex);
            GPSSignals(5,:) = xHatAux3(:,JNRIndex,bandwidthIndex);
            for i = 1:5
                sspace = 0;                 % Search space were the results will be stored
                for ii = 1:K
                    y = GPSSignals(i, (ii - 1) * Nc + (1:Nc) );   % use just 1 period of code at the time
                    % Compute the search space for a single coherent integration epoch
                    [Tsspace, searchSpace] = DftParallelCodePhaseAcquisition( y, locC, Nc, Nd, DopStep, fs, fi );
                    sspace = sspace + Tsspace;  % Non-coherently accumulate the results
                end
                
                [~, DopInd] = max(max(sspace.'));
                [~, codInd] = max(max(sspace));
                
                corrOut(i,JNRIndex,bandwidthIndex,loopIndex) = searchSpace(DopInd,codInd);
                
                correctPeak = sspace(idxFreq,idxDelay);
                                
                sspaceAux = sspace(:);
                sspaceAux = sspaceAux(sspaceAux~=correctPeak);
                
                peakAux = (mean(correctPeak - sspaceAux)^2);
                varSspace = var(sspace(:));
                generalisedSNR(i,JNRIndex,bandwidthIndex,loopIndex) =  peakAux/varSspace;
               
            end
        end
    end
end

averageGenSNR = mean(generalisedSNR, 4);

M = 20;
for bandwidthIndex = 1:length(bandwidthVector)
    for JNRIndex = 1:length(JNRVector)
        for i = 1:5
            corrBuffer = buffer(corrOut(i,JNRIndex,bandwidthIndex,:), M, 0, 'nodelay');
            Pn = sum(real(corrBuffer), 1).^2 + sum(imag(corrBuffer), 1).^2;
            Pw = sum(real(corrBuffer).^2 + imag(corrBuffer).^2, 1);
            s(i,JNRIndex,bandwidthIndex) = mean(Pn ./ Pw);
            stdS(i,JNRIndex,bandwidthIndex) = std(Pn ./ Pw);
            if s(i,JNRIndex,bandwidthIndex) < 1
                s(i,JNRIndex,bandwidthIndex) = 1;
            end
            C_N0(i,JNRIndex,bandwidthIndex) = 10*log10((1/1e-3) * ((s(i,JNRIndex,bandwidthIndex) - 1))/((M - s(i,JNRIndex,bandwidthIndex))));
        end
    end
end


for bandwidthIndex = 1:length(bandwidthVector)
    figure;
    plot(JNRVector, pow2db(averageGenSNR(1,:,bandwidthIndex)), 'Color', [0 0.4470 0.7410]);
    hold on;
    plot(JNRVector, pow2db(averageGenSNR(2,:,bandwidthIndex)), 'Color', [0.850  0 0.3250 0.0980]);
    plot(JNRVector, pow2db(averageGenSNR(3,:,bandwidthIndex)), 'Color', [0.9290 0.6940 0.1250]);
    plot(JNRVector, pow2db(averageGenSNR(4,:,bandwidthIndex)), 'Color', [0.4940 0.1840 0.5560]);
    plot(JNRVector, pow2db(averageGenSNR(5,:,bandwidthIndex)), 'Color', [0.4660 0.6740 0.1880]);
    grid on;
    ylim([-10 50]);
    ylabel('Generalised SNR [dB]')
    xlabel('JNR [dB]');
    legend('NMF-stft', 'Kalman-Perf', 'Kalman-Est.', 'NMF-fsst (Kaiser)', 'NMF-stft (Kaiser)', 'location', 'best');
%     formatFig(gcf, [dataPath 'av_snr_' num2str(bandwidthIndex)], 'en', figProp);
end

for bandwidthIndex = 1:length(bandwidthVector)
    figure;
    plot(JNRVector, squeeze(C_N0(1,:,bandwidthIndex)), 'Color', [0 0.4470 0.7410]);
    hold on;
    plot(JNRVector, squeeze(C_N0(2,:,bandwidthIndex)), 'Color', [0.8500 0.3250 0.0980]);
    plot(JNRVector, squeeze(C_N0(3,:,bandwidthIndex)), 'Color', [0.9290 0.6940 0.1250]);
    plot(JNRVector, squeeze(C_N0(4,:,bandwidthIndex)), 'Color', [0.4940 0.1840 0.5560]);
    plot(JNRVector, squeeze(C_N0(5,:,bandwidthIndex)), 'Color', [0.4660 0.6740 0.1880]);
    grid on;
    ylabel('$C/N_0$ [dB]')
    xlabel('JNR [dB]');
    legend('NMF-stft', 'Kalman-Perf', 'Kalman-Est.', 'NMF-fsst (Kaiser)', 'NMF-stft (Kaiser)', 'location', 'best');
end


rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'Misc'])
rmpath(['..' filesep 'signalsGeneration' filesep]);