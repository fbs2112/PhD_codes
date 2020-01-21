clear;
clc;
close all;

%------System parameters---------------------------------------------------
paramsSignal.Freqsamp = 32.768e6;
paramsSignal.Freqmedi = 0;
paramsSignal.FreqDopp = 1e3;
paramsSignal.Freqcode = 1.023e6;
paramsSignal.Intetime = 0.1e-3;
paramsSignal.Intenumb = fix(paramsSignal.Freqsamp*paramsSignal.Intetime);
paramsSignal.Freqwordmedi = paramsSignal.Freqmedi*2^32/paramsSignal.Freqsamp;
paramsSignal.FreqwordDopp = paramsSignal.FreqDopp*2^32/paramsSignal.Freqsamp;
paramsSignal.Freqwordcode = paramsSignal.Freqcode*2^32/paramsSignal.Freqsamp;
paramsSignal.Segnumb = 8;
paramsSignal.Nonesegment = paramsSignal.Intenumb/paramsSignal.Segnumb;
paramsSignal.Loopnumb = 50;
%--------------------------------------------------------------------------

paramsSignal.numberOfGPSSignals = 5;
paramsSignal.codeLength = 1023;

save(['.' filesep 'sim_params' filesep 'sim_params_' '3' '.mat'], 'paramsSignal'); 