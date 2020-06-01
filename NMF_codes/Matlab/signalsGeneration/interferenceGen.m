function [interferenceSignal, iflaw]  = interferenceGen(params)

Endtimelastloop = -1;
Timeofthisloop = Endtimelastloop + (1:params.Intenumb);        
% Interference signal
Phasemod = mod(Timeofthisloop, params.Noneperiod)+1;
iflaw = (params.foneperiod(Phasemod)/params.Freqsamp)';
[interferenceSignal, ~] = fmodany(iflaw);
