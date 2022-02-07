function [GPSSignals, Codesample] = GPSGen(params, varargin)
code = codeGen(params.numberOfGPSSignals, params.codeLength);

Freqwordcode = params.Freqcode*2^32/params.Freqsamp;
Endtimelastloop = -1;
Codeadderend = -1*Freqwordcode*ones(params.numberOfGPSSignals,1);
Timeofthisloop = Endtimelastloop + (1:params.Intenumb);               
Carrphase = mod(2*pi*(params.Freqmedi+params.FreqDopp)*Timeofthisloop/params.Freqsamp,2*pi);
Carrier = exp(1i*Carrphase);

% Code samples
Codeadder = mod(Codeadderend*ones(1, params.Intenumb) + (ones(params.numberOfGPSSignals, 1)*Freqwordcode)*(1:params.Intenumb), params.codeLength*2^32);

Codeadderinte = fix(Codeadder/2^32)+1;
Codesample = zeros(params.numberOfGPSSignals, params.Intenumb);

for Chanindex = 1:params.numberOfGPSSignals
    Codesample(Chanindex,:) = code(Chanindex, Codeadderinte(Chanindex,:)); 
end

if nargin == 2 && strcmp(varargin{1}, 'hilbert')
    Codesample = hilbert(Codesample);
end

% GPS signals
GPSSignals = (Codesample.*(ones(params.numberOfGPSSignals, 1)*Carrier)).';
