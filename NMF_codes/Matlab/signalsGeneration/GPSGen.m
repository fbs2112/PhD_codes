function  GPSSignals = GPSGen(params)

code = codeGen(params.numberOfGPSSignals, params.codeLength);

Endtimelastloop = -1;
Codeadderend = -1*params.Freqwordcode*ones(params.numberOfGPSSignals,1);
Timeofthisloop = Endtimelastloop + (1:params.Intenumb);               
Carrphase = mod(2*pi*(params.Freqmedi+params.FreqDopp)*Timeofthisloop/params.Freqsamp,2*pi);
Carrier = exp(1i*Carrphase);

% Code samples
Codeadder = mod(Codeadderend*ones(1, params.Intenumb) + (ones(params.numberOfGPSSignals, 1)*params.Freqwordcode)*(1:params.Intenumb), params.codeLength*2^32);
Codeadderinte = fix(Codeadder/2^32)+1;
Codesample = zeros(params.numberOfGPSSignals, params.Intenumb);

for Chanindex = 1:params.numberOfGPSSignals
    Codesample(Chanindex,:) = code(Chanindex, Codeadderinte(Chanindex,:)); 
end

% GPS signals
GPSSignals = (Codesample.*(ones(params.numberOfGPSSignals, 1)*Carrier)).';

