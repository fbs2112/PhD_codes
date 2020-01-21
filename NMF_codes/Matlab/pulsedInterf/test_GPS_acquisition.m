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
dataPath = ['..' filesep 'figs' filesep '11-12' filesep];

load sim_params_2.mat;

load(['.' filesep 'data' filesep 'nmf_testing_08.mat']);

varN = 2;
Pfa = 0.15;
thresholdBorio = sqrt(-2*varN*log(Pfa));
fi = 0;
fs = paramsSignal.Freqsamp;
Nd = 81;
DopStep = 125; %(Hz)

nbits = [2 4 8 16];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Generate the local code and resample it                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nbitsIndex = 1;
locBipolar = codeGen(1, paramsSignal.codeLength); %my GPS gen code

fc = 1.023e6;   % Code rate of the code
Tc = paramsSignal.Intetime;     % Coherent integration time in ms
Nc = floor(Tc * fs);
paramsSignal.numberOfGPSSignals = 1;

[~, locC] = GPSGen(paramsSignal);   % Resample the code at data sampling frequency
locC = quantise_gps(locC, nbits(nbitsIndex));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Main acquisition loop - evaluation of the search space         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K = 1;  % Number of non-coherent integrations
codeDl = (0:(Nc - 1)) / fs; % Code delays

JNRVector = [-15 -10 -5 10 30 50];

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
        
        signalPower = GPSSignals(i,:)*GPSSignals(i,:)'/length(GPSSignals(i,:));
        gpsSignalPowerHat = mean(abs(real(GPSSignals(i,:))))^2;
        
        CN0_SNV(JNRIndex,i) = (gpsSignalPowerHat/(signalPower - gpsSignalPowerHat))/Tc;
        M4 = mean(abs(GPSSignals(i,:)).^4);
        CN0_MM(JNRIndex,i) = (sqrt(2*(signalPower^2) - M4)/(signalPower - gpsSignalPowerHat))/Tc;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                               Plot results                                     %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        figure
        surf( codeDl*1e6, ((1:Nd) - ceil(Nd/2))*DopStep, sspace./max(sspace(:)), 'EdgeColor', 'none');
        axis tight
        xlabel('Delay [$\mu$s]')
        ylabel('Doppler Frequency (Hz)')
        zlabel('Normalized Correlation')
        if i == 1
            title(['NMF ' num2str(JNRVector(JNRIndex)) ' dB']);
        elseif i == 2
            title(['Pulse Blanking ' num2str(JNRVector(JNRIndex)) ' dB']);
        else
            title(['No Mitigation Technique ' num2str(JNRVector(JNRIndex)) ' dB']);
        end
        
        figProp = struct('size' , fontsize , 'font' , fontname , 'lineWidth' , linewidth, 'figDim', [1 1 1280 720]);
        % formatFig(gcf, [dataPath  'acq_corr'], 'en', figProp);
        
        figProp = struct('size' , fontsize , 'font' , fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
        
        % Find the Doppler frequency
        if bitand( Nd, 1 ) == 1,    % It is an odd number
            Freq = (-((Nd - 1) / 2):((Nd - 1) / 2) ) * DopStep;
        else
            Freq = (-(Nd/2):( (Nd-2) / 2 ) ) * DopStep;
        end
                        
        [maxVal, DopInd] = max(max(sspace.'));
        peakRatio(JNRIndex,i) = 2*length(GPSSignals(i,:))*maxVal/signalPower;
        [peakVector, idx] = sort(sspace(DopInd,:), 'descend');
        idxAux = idx(2:end);
        minPeakDistance = round(fs*1e-6);
        secondPeakIndexAux = find(abs(idxAux - idx(1)) > minPeakDistance, 1, 'first');
        peakRationBorre(JNRIndex,i) = maxVal/sspace(DopInd,idxAux(secondPeakIndexAux));

        DopFreq(JNRIndex, i) = Freq(DopInd);
%         figure;
%         plot( codeDl*1e6, sspace(DopInd, :) / maxVal )
%         hold on;
%         secondPeakAxis = zeros(1, length(codeDl));
%         secondPeakAxis(idxAux(secondPeakIndexAux)) = sspace(DopInd, idxAux(secondPeakIndexAux)) / maxVal;
%         plot( codeDl*1e6, secondPeakAxis, 'r' )
%         axis tight
%         grid on
%         xlabel('Delay [$\mu$s]')
%         ylabel('Normalized Correlation')
        % formatFig(gcf, [dataPath  'acq_corr_delay'], 'en', figProp);
        
        
        % Find the code delay (in chips)
        
        [maxVal, codInd] = max(max(sspace));
        delay(JNRIndex, i) = codeDl(codInd)*1e6;
        codeDelay(JNRIndex, i) = 1023 - (codInd - 1) / fs * fc;
        % codeDelay = codInd / fs * fc
%         figure
%         plot(Freq, sspace(:, codInd) / maxVal)
%         axis tight
%         grid on
%         xlabel('Doppler Frequency (Hz)')
%         ylabel('Normalized Correlation')
        % formatFig(gcf, [dataPath  'acq_corr_doppler'], 'en', figProp);
    end
end

sspace = 0;

signalPower = locC*locC'/length(locC); 
gpsSignalPowerHat = mean(abs(real(signalPower)))^2;
CN0_SNV(JNRIndex + 1,i + 1) = (gpsSignalPowerHat/(signalPower - gpsSignalPowerHat))/Tc;

M4 = mean(abs(GPSSignals(i,:)).^4);
CN0_MM(JNRIndex+1,i+1) = (sqrt(2*signalPower^2 - M4)/(signalPower - gpsSignalPowerHat))/Tc;
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
title(['GPS Signal ' num2str(JNRVector(JNRIndex)) ' dB']);

figure;
plot(JNRVector, zeroedSamples)
xlabel('JNR [dB]')
ylabel('Discarded samples');

rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'Misc'])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);