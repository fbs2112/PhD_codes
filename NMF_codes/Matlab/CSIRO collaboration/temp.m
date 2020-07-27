clear;
clc;
close all;

for i = 1:198
    load(['.' filesep 'data' filesep 'dataParks_' num2str(i) '.mat']);
    parkesSignal = parksSignalAux;
    save(['.' filesep 'data' filesep 'dataParkes_' num2str(i) '.mat'], 'parkesSignal', 'Fs');
end