clear;
clc;
close all;

fileName = 'fs12mhz_2prn_1172fc_144s';
matObj = matfile(['.' filesep 'IRNSS signals' filesep fileName '.mat']);

fs = 12e6;
IF = 4.092e6;

numberOfLoops = 48/100e-3;
lengthData = floor(fs*100e-3);
t = 0:1/fs:100e-3 -1/fs;
carrierConj = exp(-1j*2*pi*(IF)*t).';
lp_filter = dsp.LowpassFilter('SampleRate', fs, 'PassbandFrequency', 2e6, 'StopbandFrequency', 2.5e6);

for i = 1:numberOfLoops
    i
   
    samplesAux = (matObj.samples( (i-1)*lengthData + 1:i*lengthData,1));
    samplesAux = samplesAux.*carrierConj;
    samplesAux = lp_filter(samplesAux);
    samples(:,i) = samplesAux;
    
end
samples = samples(:);

save(['.' filesep 'IRNSS signals' filesep fileName '_downconverted_low_pass.mat'], 'samples', 'fs', '-v7.3');