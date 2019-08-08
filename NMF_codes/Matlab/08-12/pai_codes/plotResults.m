clear;
clc;
close all;

load resultsPai01.mat

signalIndex = 1:50;
noSignalIndex = 60:273;

x = mean(GoFBlockPd(signalIndex));
y = mean(GoFBlockPd(noSignalIndex));