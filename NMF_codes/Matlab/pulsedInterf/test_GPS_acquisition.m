clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');

addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
addpath(['..' filesep 'Misc'])

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
dataPath = ['..' filesep 'figs' filesep '11-12' filesep];

load sim_params_2.mat;

load(['.' filesep 'data' filesep 'nmf_testing_03.mat']);

JNRVector = [-10 0 10 20 30];
JNRIndex = 1;
varN = 2;
Pfa = 0.15;
thresholdBorio = sqrt(-2*varN*log(Pfa));
xHatBorio = mixtureSignal(:,JNRIndex);
xHatBorio(abs(xHatBorio) > thresholdBorio) = 0;

GPSSignals(1,:) = xHat(:,2, JNRIndex).';
GPSSignals(2,:) = xHatBorio.';
GPSSignals(3,:) = mixtureSignal(:,JNRIndex).';
fi = 0;
fs = paramsSignal.Freqsamp; 
Nd = 81; 
DopStep = 50; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Generate the local code and resample it                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sat = 7;       % PRN of the satellite to acquire
loc = codeGen(1, paramsSignal.codeLength); %my GPS gen code

% loc = GpsCaCode( sat ); % The code is in binary format - need to be converted
                        % into a bipolar format
locBipolar = 1 - 2 * loc;   % Code in bipolar format

fc = 1.023e6;   % Code rate of the code
Tc = paramsSignal.Intetime;     % Coherent integration time in ms
Nc = floor(Tc * fs);
locC = ResampleCode( locBipolar, Nc, fs, 0, fc );   % Resample the code at
                                                    % data sampling
                                                    % frequency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Main acquisition loop - evaluation of the search space         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K = 1;  % Number of non-coherent integrations
sspace = 0;                 % Search space were the results will be stored
codeDl = (0:(Nc - 1)) / fs; % Code delays

for i = 1:3
    for ii = 1:K
        y =  	GPSSignals(i, (ii - 1) * Nc + (1:Nc) );   % use just 1 period of code at the time    
        % Compute the search space for a single coherent integration epoch
        Tsspace = DftParallelCodePhaseAcquisition( y, locC, Nc, Nd, DopStep, fs, fi );

        sspace = sspace + Tsspace;  % Non-coherently accumulate the results 
    end 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                             Add decision logic here                            % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Noise floor estimation and threshold evaluation

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                               Plot results                                     % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure
    surf( codeDl*1e6, ((1:Nd) - ceil(Nd/2))*DopStep, sspace./max(sspace(:)), 'EdgeColor', 'none');
    axis tight
    xlabel('Delay [$\mu$s]')
    ylabel('Doppler Frequency (Hz)')
    zlabel('Normalized Correlation')
    figProp = struct('size' , fontsize , 'font' , fontname , 'lineWidth' , linewidth, 'figDim', [1 1 1280 720]);
    % formatFig(gcf, [dataPath  'acq_corr'], 'en', figProp);

    figProp = struct('size' , fontsize , 'font' , fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);

    % Find the Doppler frequency
    figure
    if bitand( Nd, 1 ) == 1,    % It is an odd number 
        Freq = (-((Nd - 1) / 2):((Nd - 1) / 2) ) * DopStep;
    else
        Freq = (-(Nd/2):( (Nd-2) / 2 ) ) * DopStep;
    end

    [maxVal DopInd] = max(max(sspace.'));
    DopFreq = Freq(DopInd)
    plot( codeDl*1e6, sspace(DopInd, :) / maxVal )
    axis tight
    grid on
    xlabel('Delay [$\mu$s]')
    ylabel('Normalized Correlation')
    % formatFig(gcf, [dataPath  'acq_corr_delay'], 'en', figProp);


    % Find the code delay (in chips)
    figure;
    [maxVal codInd] = max(max(sspace));
    delay = codeDl(codInd)*1e6
    codeDelay = 1023 - (codInd - 1) / fs * fc
    % codeDelay = codInd / fs * fc

    plot(Freq, sspace(:, codInd) / maxVal)
    axis tight
    grid on
    xlabel('Doppler Frequency (Hz)')
    ylabel('Normalized Correlation')
    % formatFig(gcf, [dataPath  'acq_corr_doppler'], 'en', figProp);
end

rmpath(['..' filesep 'Misc'])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);