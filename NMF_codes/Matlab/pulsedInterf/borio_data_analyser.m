clear;
clc;
close all;

addpath(['.' filesep 'data']);
addpath(['..' filesep 'Sigtools' filesep])
addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
addpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])

load sim_params_1.mat;


load raw_data_borio.mat;
rawdata = rawdata(:);

params.fs = 10e6;
params.nfft = 256;
params.nperseg = 256;
params.overlap = params.nperseg - 1;
params.hop_size = params.nperseg - params.overlap;
params.window = ones(params.nperseg, 1);
params.specType = 'power';
params.numberOfSources = 6;
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 1000;
params.tolChange = 1e-6;
params.tolError = 1e-6;
params.repetitions = 1;
params.JNRVector = 0;

[WInterf, HInterf, error, PxxInterf, f, t] = nmf_eval_v2(rawdata, params);

figure;
surf(t*1e6, f/1e6, 10*log10(abs(PxxInterf{1,1}).^2), 'EdgeColor', 'none');
axis xy;
view(0, 90);
xlabel('Time [$\mu$s]');
ylabel('Frequency [MHz]');
xlim([0 max(t)]*1e6);
ylim([min(f) max(f)]/1e6);
ax = gca;
% set(ax, 'CLim', [-20 40], 'colormap', jet);
set(ax, 'colormap', jet);

c = colorbar;
c.Label.String = '[dB]';
c.TickLabelInterpreter = 'latex';

figure;
plot(error);

figure;
plot(f/1e6, WInterf{1,1})

figure;
plot(t*1e6, HInterf{1,1})

rmpath(['.' filesep 'data']);addpath(['.' filesep 'data']);
rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
rmpath(['..' filesep 'Sigtools' filesep 'NMF_algorithms'])
