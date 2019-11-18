clear;
clc;
close all;


addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);

load sim_params_2.mat;


load(['.' filesep 'data' filesep 'data_NMF_3.mat']);
load(['.' filesep 'data' filesep 'data_no_processing_3.mat']);

GPSSignals = xHat(:,2).';
% GPSSignals = mixtureSignal.';

% GPSSignals = GPSGen(paramsSignal).';

% GPSSignals = GPSSignals(1:numberOfRawSamples,1);
% GPSSignalsPower = pow_eval(GPSSignals);
fi = 0;
fs = paramsSignal.Freqsamp; 
Nd = 81; 
DopStep = 125; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Generate the local code and resample it                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sat = 7;       % PRN of the satellite to acquire
loc = codeGen(1, paramsSignal.codeLength); %my GPS gen code

% loc = GpsCaCode( sat ); % The code is in binary format - need to be converted
                        % into a bipolar format
locBipolar = 1 - 2 * loc;   % Code in bipolar format

fc = 1.023e6;   % Code rate of the code
Tc = 0.0001;     % Coherent integration time in ms
Nc = floor(Tc * fs);
locC = ResampleCode( locBipolar, Nc, fs, 0, fc );   % Resample the code at
                                                    % data sampling
                                                    % frequency
locC = locC(1:length(GPSSignals));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Main acquisition loop - evaluation of the search space         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K = 1;  % Number of non-coherent integrations
sspace = 0;                 % Search space were the results will be stored
codeDl = (0:(Nc - 1)) / fs; % Code delays

for ii = 1:K,
    y =  	GPSSignals( (ii - 1) * Nc + (1:Nc) );   % use just 1 period of code at the time
    
%     y =  	GPSSignals;   % use just 1 period of code at the time
    
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
surf( codeDl, ((1:Nd) - ceil(Nd/2))*DopStep, sspace, 'EdgeColor', 'none');
axis tight
set( gca, 'FontSize', 16 )
xlabel('Code delay [samples]')
ylabel('Frequency [Hz]')

% Find the Doppler frequency
figure
if bitand( Nd, 1 ) == 1,    % It is an odd number 
    Freq = (-((Nd - 1) / 2):((Nd - 1) / 2) ) * DopStep;
else
    Freq = (-(Nd/2):( (Nd-2) / 2 ) ) * DopStep;
end

[maxVal DopInd] = max(max(sspace.'));
DopFreq = Freq(DopInd)
plot( codeDl, sspace(DopInd, :) / maxVal )
axis tight
grid on
set( gca, 'FontSize', 16)
xlabel('Code delay (samples)')
ylabel('Normalized Correlation')

% Find the code delay (in chips)
figure
[maxVal codInd] = max(max(sspace));
codeDelay = 1023 - (codInd - 1) / fs * fc
% codeDelay = codInd / fs * fc

plot(Freq, sspace(:, codInd) / maxVal)
axis tight
grid on
set( gca, 'FontSize', 16)
xlabel('Doppler frequency (Hz)')
ylabel('Normalized Correlation')

rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);