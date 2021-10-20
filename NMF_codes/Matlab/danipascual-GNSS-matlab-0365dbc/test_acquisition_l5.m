clear;
clc;
close all;


addpath(['.' filesep 'source']);
addpath(['.' filesep 'source' filesep 'functions' filesep 'signal_processing']);
addpath(['.' filesep 'prn_codes']);
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
addpath(['..' filesep 'Misc'])
addpath(['..' filesep 'Sigtools' filesep])

load sim_params_3.mat;

fi = 0;
fs = 50e6;
Nd = 81;
DopStep = 125; %(Hz)

delay = 10e-6;
trueDoppler = 1.2e3;
paramsSignal.numberOfGPSSignals = 1;

[GPSL5_I, GPSL5_Q] = GNSSsignalgen(1, 'L5', fs, 1);   % Resample the code at data sampling frequency
% GPSL5_I = [GPSL5_I(end - round(delay*fs)+1:end,:);GPSL5_I(1:end - round(delay*fs),:)]; % Introducing artificial code delay
GPSL5 = GPSL5_I + 1j*GPSL5_Q;
GPSL5 = GPSL5.';

Tc = paramsSignal.Intetime;     % Coherent integration time in ms
Nc = floor(Tc * fs);

K = 1;  % Number of non-coherent integrations
codeDl = (0:(Nc - 1)) / fs; % Code delays

if bitand( Nd, 1 ) == 1    % It is an odd number
    Freq = (-((Nd - 1) / 2):((Nd - 1) / 2) ) * DopStep;
else
    Freq = (-(Nd/2):( (Nd-2) / 2 ) ) * DopStep;
end

Timeofthisloop = 0:(fs*1e-3)-1;

Carrphase = mod(2*pi*(1e3)*Timeofthisloop/fs,2*pi);
Carrier = exp(1i*Carrphase);

for i = 1:1
    sspace = 0;                 % Search space were the results will be stored
    for ii = 1:K
        [I, Q] = GNSSsignalgen(1, 'L5', fs, 1);    % use just 1 period of code at the time
        y = (I + 1j*Q).';
        y = y.*Carrier;
        y = [y(end - round(delay*fs)+1:end) y(1:end - round(delay*fs))]; % Introducing artificial code delay

        % Compute the search space for a single coherent integration epoch
        Tsspace = DftParallelCodePhaseAcquisition( y, GPSL5, Nc, Nd, DopStep, fs, fi );
        sspace = sspace + Tsspace;  % Non-coherently accumulate the results
    end
    
    figure
    surf( codeDl*1e6, ((1:Nd) - ceil(Nd/2))*DopStep, sspace./max(sspace(:)), 'EdgeColor', 'none');
    axis tight
    xlabel('Delay [$\mu$s]')
    ylabel('Doppler Frequency (Hz)')
    zlabel('Normalised Correlation')
  
    [maxVal, DopInd] = max(max(sspace.'));
    DopFreq = Freq(DopInd);
    [maxVal, codInd] = max(max(sspace));
    delayEst = codeDl(codInd)*1e6;
end

rmpath(['.' filesep 'source']);
rmpath(['.' filesep 'source' filesep 'functions' filesep 'signal_processing']);
rmpath(['.' filesep 'prn_codes']);
rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'Misc'])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);