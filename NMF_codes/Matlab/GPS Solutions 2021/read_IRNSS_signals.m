clear;
clc;
close all;

matObj = matfile(['.' filesep 'IRNSS signals' filesep 'fs12mhz_2prn_1172fc_48s.mat']);

fs = 12e6;
samples = matObj.samples(1:fs*100e-3, 1);
length_data = length(samples)/fs;

PRN_ID = 2;

IRNSS_signals = IRNSS_PRN(PRN_ID, 'L5');

fc = 1.023e6;   % Code rate of the code
Tc = 0.001;     % Coherent integration time in ms
Nc = Tc * fs;
IRNSS_signals_resampled = ResampleCode(IRNSS_signals, Nc, fs, 0, fc); 

figure;
pspectrum(IRNSS_signals_resampled, fs);

figure; 
pspectrum(samples, fs);

figure;
[pxx,f] = pwelch(samples,[],[],[],fs,'centered', 'power');
plot(f,10*log10(abs(pxx)));
% semilogy(f,10*log10(abs(p)))
IF = 4.092e6;

t = 0:1/fs:length_data-1/fs;
% Carrphase = mod(2*pi*(IF)*t,2*pi);
carrierConj = exp(-1j*2*pi*(IF)*t).';

downConvertedSignal = samples.*carrierConj;
figure;
[pxx,f] = pwelch(downConvertedSignal,[],[],[],fs,'centered');
plot(f,10*log10(abs(pxx)));


lp_filter = dsp.LowpassFilter('SampleRate', fs, 'PassbandFrequency', 2e6, 'StopbandFrequency', 2.5e6);

filteredSignal  = lp_filter(downConvertedSignal);
% fvtool(lp_filter, 'Analysis', 'Freq');

figure;
[pxx,f] = pwelch(filteredSignal,[],[],[],fs,'centered');
plot(f,10*log10(abs(pxx)));

