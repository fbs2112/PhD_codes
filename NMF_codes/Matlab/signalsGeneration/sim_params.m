clear;
clc;
close all;

%------System parameters---------------------------------------------------
paramsSignal.Freqsamp = 32.768e6;
paramsSignal.Freqmedi = 0;
paramsSignal.FreqDopp = 1e3;
paramsSignal.Freqcode = 1.023e6;
paramsSignal.Intetime = 1e-3;
paramsSignal.Intenumb = fix(paramsSignal.Freqsamp*paramsSignal.Intetime);
paramsSignal.Freqwordmedi = paramsSignal.Freqmedi*2^32/paramsSignal.Freqsamp;
paramsSignal.FreqwordDopp = paramsSignal.FreqDopp*2^32/paramsSignal.Freqsamp;
paramsSignal.Freqwordcode = paramsSignal.Freqcode*2^32/paramsSignal.Freqsamp;
paramsSignal.Segnumb = 8;
paramsSignal.Nonesegment = paramsSignal.Intenumb/paramsSignal.Segnumb;
paramsSignal.Loopnumb = 50;
%--------------------------------------------------------------------------



%------Interference parameters-------------------------------------------
% paramsSignal.Sweeptimeup = 6.8e-6;                             % ramp-up sweep time
% paramsSignal.Sweeptimedown = 1.9e-6;                           % ramp-down sweep time
% paramsSignal.Sweeptime = 8.7e-6;                               % sweep time
% paramsSignal.Nup = round(paramsSignal.Sweeptimeup*paramsSignal.Freqsamp);               % number of samples with a ramp-up sweep time 
% paramsSignal.Ndown = round(paramsSignal.Sweeptimedown*paramsSignal.Freqsamp);           % number of samples with a ramp-down sweep time
% paramsSignal.Noneperiod = paramsSignal.Nup+paramsSignal.Ndown;                         % number of samples with a sweep time
% paramsSignal.IFmin = 4e6;                                      % start frequency
% paramsSignal.IFmax = 8e6;                                      % end frequency
% 
% paramsSignal.foneperiod(1:paramsSignal.Noneperiod) = linspace(paramsSignal.IFmin, paramsSignal.IFmax, paramsSignal.Noneperiod);
% paramsSignal.Initphase = 0;
%-------------------------------------------------------------------------


paramsSignal.numberOfGPSSignals = 5;
paramsSignal.codeLength = 1023;

save(['.' filesep 'sim_params' filesep 'sim_params_' '1' '.mat'], 'paramsSignal'); 