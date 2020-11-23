clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['.' filesep 'data']);
addpath(['..' filesep '.' filesep 'Misc'])

linewidth = 1.5;
fontname = 'Times';
fontsize = 30;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
dataPath = ['..' filesep '.' filesep 'figs' filesep '2020-11-05' filesep];

load ADSB_label.mat;
load template_ADSB.mat;

params.fs = 128e6;
params.nfft = 256;
params.nperseg = 256;
params.overlap = params.nperseg/2;
params.hop_size = params.nperseg - params.overlap;
params.window = ones(params.nperseg, 1);
params.specType = 'power';

fVector = linspace(0, params.fs - params.fs/params.nfft, params.nfft);
centralFrequency = 66e6;

template = zeros(length(fVector), 1);
template(128:138) = 1;
figure;
plot(fVector/params.fs * 2, template);
xlabel('Normalized Frequency [$\times \pi$ rad/sample]');
ylabel('Magnitude');
% xlim([min(f) max(f)]/1e6 + 1152);
ylim([0 1.1]);
formatFig(gcf, [dataPath 'template'], 'en', figProp);



rmpath(['..' filesep '.' filesep 'Misc'])
rmpath(['.' filesep 'data']);
rmpath(['..' filesep '.' filesep 'Sigtools' filesep])
rmpath(['..' filesep  '.' filesep 'Sigtools' filesep 'NMF_algorithms'])