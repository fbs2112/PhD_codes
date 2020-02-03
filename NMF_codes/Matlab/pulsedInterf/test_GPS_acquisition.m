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

load(['.' filesep 'data' filesep 'nmf_testing_11.mat']);

varN = 2;
Pfa = 0.15;
thresholdBorio = sqrt(-2*varN*log(Pfa));
fi = 0;
fs = paramsSignal.Freqsamp;
Nd = 81;
DopStep = 125; %(Hz)

trueDelay = 10e-6;
trueDoppler = 1.5e3;
paramsSignal.numberOfGPSSignals = 1;

[~, locC] = GPSGen(paramsSignal);   % Resample the code at data sampling frequency

JNRVector = [-15 -10 -5 10 30 50];

fc = 1.023e6;   % Code rate of the code
Tc = paramsSignal.Intetime;     % Coherent integration time in ms
Nc = floor(Tc * fs);

K = 1;  % Number of non-coherent integrations
codeDl = (0:(Nc - 1)) / fs; % Code delays

nbits = [2 4 8 16];

for nbitsIndex = 3:3%length(nbits)
    locC = quantise_gps(locC, nbits(nbitsIndex));
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
            peakRatio(JNRIndex,i) = maxVal/signalPower;
            [peakVector, idx] = sort(sspace(DopInd,:), 'descend');
            idxAux = idx(2:end);
            minPeakDistance = round(fs*1e-6);
            secondPeakIndexAux = find(abs(idxAux - idx(1)) > minPeakDistance, 1, 'first');
            peakRatioBorre(JNRIndex,i) = maxVal/sspace(DopInd,idxAux(secondPeakIndexAux));
            
            DopFreq(JNRIndex,i) = Freq(DopInd);
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
            
            
            
            %         Carrphase = mod(2*pi*(trueDoppler)*codeDl,2*pi);
            %         Carrier = exp(1i*Carrphase).';
            %         GPSSignalNoDoppler = (GPSSignals(i,:).').*conj(Carrier);
            %         locCDelayed = [locC(end - round(trueDelay*fs)+1:end) locC(1:end - round(trueDelay*fs))]; % Introducing artificial code delay
            %         navSignalHat = GPSSignalNoDoppler.*(locCDelayed).';
            %
            %         signalPower = navSignalHat'*navSignalHat/length(navSignalHat);
            %         gpsSignalPowerHat = mean(abs(real(navSignalHat)))^2;
            %         CN0_SNV(JNRIndex,i) = pow2db((gpsSignalPowerHat/(signalPower - gpsSignalPowerHat))/Tc);
            %
            %         Carrphase = mod(2*pi*(DopFreq(JNRIndex,i))*codeDl,2*pi);
            %         Carrier = exp(1i*Carrphase).';
            %         GPSSignalNoDoppler = (GPSSignals(i,:).').*conj(Carrier);
            %         locCDelayed = [locC(end - round(delay(JNRIndex,i)*fs/1e6)+1:end) locC(1:end - round(delay(JNRIndex,i)*fs/1e6))]; % Introducing artificial code delay
            %         navSignalHat = GPSSignalNoDoppler.*(locCDelayed).';
            %
            %         signalPower = navSignalHat'*navSignalHat/length(navSignalHat);
            %         gpsSignalPowerHat = mean(abs(real(navSignalHat)))^2;
            %         CN0_SNV2(JNRIndex,i) = pow2db((gpsSignalPowerHat/(signalPower - gpsSignalPowerHat))/Tc);
        end
    end
    
end
% sspace = 0;
%
% signalPower = locC*locC'/length(locC);
% gpsSignalPowerHat = mean(abs(real(signalPower)))^2;
% CN0_SNV(JNRIndex + 1,i + 1) = pow2db(gpsSignalPowerHat/(signalPower - gpsSignalPowerHat))/Tc;

figure;
plot(JNRVector, zeroedSamples)
xlabel('JNR [dB]')
ylabel('Discarded samples');

rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'Misc'])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);