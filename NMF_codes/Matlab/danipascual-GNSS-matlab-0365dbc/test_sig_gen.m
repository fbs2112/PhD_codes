clear;
clc;
close all;

addpath(['.' filesep 'source']);
addpath(['.' filesep 'source' filesep 'functions' filesep 'signal_processing']);
addpath(['.' filesep 'prn_codes']);
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);

load sim_params_3.mat;

fs = 32.768e6;

fs = 50e6;

[I, Q] = GNSSsignalgen(1, 'L5', fs, 1);
[I2, Q2] = GNSSsignalgen(1, 'L5I', fs, 1);
% [~,GPSSignals] = GPSGen(paramsSignal);

figure;
periodogram(I + 1j*Q,rectwin(length(I2)),length(I2),fs)
figure;
periodogram(I2,rectwin(length(I2)),length(I2),fs)


rmpath(['.' filesep 'source']);
rmpath(['.' filesep 'source' filesep 'functions' filesep 'signal_processing']);
rmpath(['.' filesep 'prn_codes']);
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);