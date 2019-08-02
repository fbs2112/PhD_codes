clear;
clc;
close all

addpath(['..' filesep 'Sigtools' filesep])

Freqsamp = 32.768e6;
Freqmedi = 0;
FreqDopp = 1e3;
Freqcode = 1.023e6;
Intetime = 1e-3;
Intenumb = fix(Freqsamp*Intetime);
Freqwordmedi = Freqmedi*2^32/Freqsamp;
FreqwordDopp = FreqDopp*2^32/Freqsamp;
Freqwordcode = Freqcode*2^32/Freqsamp;
Segnumb = 8;
Nonesegment = Intenumb/Segnumb;
Loopnumb = 50;


Channumb = 5;
SNR = -25*ones(Channumb,1);
JNR = 0; 


Sweeptimeup1 = 6.8e-6;                             % ramp-up sweep time
Sweeptimedown1 = 1.9e-6;                           % ramp-down sweep time
Sweeptime1 = 8.7e-6;                               % sweep time
Nup1 = round(Sweeptimeup1*Freqsamp);               % number of samples with a ramp-up sweep time 
Ndown1 = round(Sweeptimedown1*Freqsamp);           % number of samples with a ramp-down sweep time
Noneperiod1 = Nup1+Ndown1;                         % number of samples with a sweep time
IFmin1 = 4e6;                                      % start frequency
IFmax1 = 8e6;                                      % end frequency

foneperiod(1:Noneperiod1) = linspace(IFmin1,IFmax1,Noneperiod1);
Initphase = 0;

Sweeptime = Sweeptime1;
Noneperiod = Noneperiod1;




Codelen = 1023; 

% Generate PRN code
% Generate unipolar Gold code
Goldcode = zeros(Channumb,1023);                   % unipolar code
G1 = ones(Channumb,10);                            % register 1
G2 = ones(Channumb,10);                            % register 2
Goldphase = [2,6;3,7;4,8;1,8;2,9];                 % Gold code phase 
for index = 1:Codelen 
    for channelindex = 1:Channumb
        Goldcode(channelindex,index) = mod(G1(channelindex,10)+mod(sum(G2(channelindex,Goldphase(channelindex,1))+G2(channelindex,Goldphase(channelindex,2))),2),2);
    end
    % feedback
    G1feedback = [G1(:,3) G1(:,10)];
    G2feedback = [G2(:,2) G2(:,3) G2(:,6) G2(:,8) G2(:,9) G2(:,10)];
    % shift
    G1(:,[2:10]) = G1(:,[1:9]);
    G2(:,[2:10]) = G2(:,[1:9]);
    G1(:,1) = mod(sum(G1feedback,2),2);
    G2(:,1) = mod(sum(G2feedback,2),2);
end

code = 2*Goldcode-1;

Endtimelastloop = -1;
Codeadderend = -1*Freqwordcode*ones(Channumb,1);


Timeofthisloop = Endtimelastloop + [1:Intenumb];               
Carrphase = mod(2*pi*(Freqmedi+FreqDopp)*Timeofthisloop/Freqsamp,2*pi);
Carrier = exp(1i*Carrphase);
Endtimelastloop = Timeofthisloop(Intenumb);

% Code samples
Codeadder = mod(Codeadderend*ones(1,Intenumb) + (ones(Channumb,1)*Freqwordcode)*(1:Intenumb),Codelen*2^32);
Codeadderinte = fix(Codeadder/2^32)+1;
for Chanindex = 1:Channumb
    Codesample(Chanindex,:) = code(Chanindex,Codeadderinte(Chanindex,:)); 
end
Codeadderend = Codeadder(:,Intenumb);

% Desired signal
Desiredsignal = Codesample.*(ones(Channumb,1)*Carrier);
Exs = mean(abs(Desiredsignal).^2,2);

% GPS signal corrupted by additive Gaussian white noise
% Noise samples
noisearray = randn(1,Intenumb) + 1i*randn(1,Intenumb);
Exn = mean(abs(noisearray).^2,2);
hs = sqrt(Exn*10.^(SNR/10)./Exs);
Signalnointe = (sum((hs*ones(1,Intenumb)).*Desiredsignal,1) + noisearray).';

% GPS signal corrupted by interference

% Interference samples
Phasemod = mod(Timeofthisloop,Noneperiod)+1;
iflaw = (foneperiod(Phasemod)/Freqsamp)';
[int,iflawsim] = fmodany(iflaw);
Exi = mean(abs(int).^2);
hi = sqrt(Exn*10^(JNR/10)/Exi);
Signal = Signalnointe + hi*int;


% Interference-free signal within assessment window
noiseassess = randn(1,Intenumb) + 1i*randn(1,Intenumb);
pow_noise = pow_eval(noiseassess);
noiseassess = noiseassess.*sqrt(1/pow_noise);
Signalass = (sum((hs*ones(1,Intenumb)).*Desiredsignal,1) + noiseassess).';

save(['..' filesep 'data' filesep '08-12' filesep 'gnss_signals2.mat'], 'Signal', 'Signalass');

rmpath(['..' filesep 'Sigtools' filesep])
