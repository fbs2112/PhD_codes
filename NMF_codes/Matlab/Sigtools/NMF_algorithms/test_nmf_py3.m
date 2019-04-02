clear;
clc;
close all;


params.numberOfSources = 1;
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 1000;
params.random_state = 42;
params.tolChange = 1e-6;
params.tolError = 1e-6;


X = [1,2,3;4,5,6;7,8,9];
W0 = [1 1;8 2;1 3];
H0 = [1 2 2;3 3 8];


[W, H, reconstructError] = nmf_py3(X, params);

plot(20*log10(abs(reconstructError(1:find(reconstructError, 1, 'last')))));
