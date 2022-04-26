clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');

addpath(['..' filesep 'Misc'])
addpath(['..' filesep 'Sigtools']);

linewidth = 1.5;
fontname = 'Times';
fontsize = 38;
dataPath = ['..' filesep 'figs' filesep '2021-02-04' filesep];
figProp = struct('size' , fontsize , 'font' , fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);

fs = 12e6;
fi = 0;
fc = 1.023e6;   % Code rate of the code

Nd = 81;
DopStep = 125; %(Hz)

delay = 10e-6;

Tc = 1e-3;     % Coherent integration time in ms
Nc = floor(Tc * fs);

PRN_ID = 2;
IRNSS_signals = IRNSS_PRN(PRN_ID, 'L5');
IRNSS_signals_resampled = ResampleCode(IRNSS_signals, Nc, fs, 0, fc);

matObj = matfile(['.' filesep 'IRNSS signals' filesep 'fs12mhz_2prn_1172fc_48s_downconverted_low_pass.mat']);

K = 1;  % Number of non-coherent integrations
codeDl = (0:(Nc - 1)) / fs; % Code delays

if bitand( Nd, 1 ) == 1    % It is an odd number
    Freq = (-((Nd - 1) / 2):((Nd - 1) / 2) ) * DopStep;
else
    Freq = (-(Nd/2):( (Nd-2) / 2 ) ) * DopStep;
end
JNRIndex = 1;
stepSize = fs*1e-3;
% acquired_IRNSS = buffer(samples, fs*1e-3, 0, 'nodelay');
for j = 1:1%size(acquired_IRNSS, 2)
    
    for i = 1:1
        sspace = 0;                 % Search space were the results will be stored
        y = matObj.samples((j-1)*stepSize+1:j*stepSize,1);
        for ii = 1:K
            
            % Compute the search space for a single coherent integration epoch
            Tsspace = DftParallelCodePhaseAcquisition( y.', IRNSS_signals_resampled.', Nc, Nd, DopStep, fs, fi );
            sspace = sspace + Tsspace;  % Non-coherently accumulate the results
        end
        
        aux = ((1:Nd) - ceil(Nd/2))*DopStep;
        figure
        surf( codeDl*1e6, ((1:Nd) - ceil(Nd/2))*DopStep, sspace./max(sspace(:)), 'EdgeColor', 'none');
        axis tight
        xlabel('$\tau$ [$\mu$s]')
        ylabel('$F_\mathrm{d}$ [Hz]')
        zlabel('CAF')
        [maxVal, DopInd] = max(max(sspace.'));
        DopFreq(i,JNRIndex) = Freq(DopInd);
        [maxVal, codInd] = max(max(sspace));
        delayEst(i,JNRIndex) = codeDl(codInd)*1e6;
    end
%     close all;
    
end

rmpath(['..' filesep 'Misc'])
rmpath(['..' filesep 'Sigtools']);