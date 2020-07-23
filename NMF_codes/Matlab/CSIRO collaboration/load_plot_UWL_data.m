clear;
clc;
close all;

addpath(['.' filesep 'data']);

Fs = 128e6;

fid = fopen('Parkes_S04_113_065200.vdf','r');
dat_in = fread(fid,'uint16');
fclose(fid);
dat_in = dat_in - 2^15; % offset binary

tmp = (dat_in(1:2:end) + 1i*dat_in(2:2:end));
clear dat_in;

%% Average down the noise and plot the power spectrum
NFFT = 512;
N = 10000;
xt = reshape(tmp(1:N*NFFT),NFFT,N);
ft = fft(xt);
spect = sum((ft.*conj(ft)).');

pwr_spect = 20*log10(fftshift(spect));
pwr_spect = pwr_spect - max(pwr_spect);
xax = (0:NFFT-1)/NFFT*Fs + 1088e6;

figure(1)
clf
plot(xax/1e6,pwr_spect)
% xlim([1088 1216])
line([1090 1090],[-10 0],'linestyle','--','color','k')
grid on
xlabel('Frequency (MHz)')
ylabel('Power (dB)')

rmpath(['.' filesep 'data']);

%%
%save data

dataLength = round(10e-3*Fs);
parksSignalAux = buffer(tmp, dataLength, 0, 'nodelay');

for i = 1:size(parksSignalAux, 2)
    parksSignal = parksSignalAux(:,i);
    save(['.' filesep 'data' filesep 'dataParks_' num2str(i) '.mat'], 'parksSignal', 'Fs');
end