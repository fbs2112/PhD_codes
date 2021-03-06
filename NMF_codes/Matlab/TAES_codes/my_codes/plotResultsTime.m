clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '..' filesep '.' filesep 'Misc'])
if isunix
    addpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/my_results/time_results/');
    
else
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep 'time_results' filesep])
end  

load resultsTime1.mat;
disp(['NMF detector average time: ' num2str(mean(elapsedTime)*1e3) ' +- ' num2str(var(elapsedTime)*1e3) ' ms']);

if isunix
    rmpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/my_results/time_results/');
    
else
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'my_results' ...
        filesep 'time_results' filesep])
end  


if isunix
    addpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/pai_results/time_results/');
    
else
    addpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep 'time_results' filesep])

end  

load resultsTime1.mat;
disp(['Pai''s average time (L = 3): ' num2str(mean(elapsedTime)*1e3) ' +- ' num2str(var(elapsedTime)*1e3) ' ms']);
load resultsTime2.mat;
disp(['Pai''s average time (L = 19): ' num2str(mean(elapsedTime)*1e3) ' +- ' num2str(var(elapsedTime)*1e3) ' ms']);
load resultsTime3.mat;
disp(['Pai''s average time (L = 11): ' num2str(mean(elapsedTime)*1e3) ' +- ' num2str(var(elapsedTime)*1e3) ' ms']);

if isunix
    rmpath('/home/felipe/Dropbox/Doctorate/Research/data/TAES_data/new_data/pai_results/time_results/');
    
else
    rmpath(['..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep '..' filesep 'Dropbox' filesep ...
        'Doctorate' filesep 'Research' filesep 'data' filesep 'TAES_data' filesep 'new_data' filesep 'pai_results' ...
        filesep 'time_results' filesep])
end  