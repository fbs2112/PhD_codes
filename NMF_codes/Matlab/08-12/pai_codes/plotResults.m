clear;
clc;
close all;

addpath(['.' filesep 'data']);

load resultsPai01.mat

PfaVector = [1e-5 1e-4 1e-3 1e-2 1e-1 1];

for i = 1:length(PfaVector)
   x = squeeze(GoFBlockDeteflag(i,:,:));
   tpRegion = x(28,:);
   fpRegion = x(70:150,:);
   fpRegion = fpRegion(:);
   tp(i) = sum(tpRegion);
   fn(i) = length(tpRegion) - tp(i);
   
   fp(i) = sum(fpRegion);
   tn(i) = length(fpRegion) - fp(i);
end

tpr = tp./(tp+fn);
fpr = fp./(fp+tn);

figure;
plot(fpr, tpr);
hold on;
xAux = linspace(0, 1, numel(fpr));
yAux = linspace(0, 1, numel(fpr));
plot(xAux, yAux, '--');
box 'off'
grid on
xlabel('Probability of false alarm');
ylabel('Probability of detection');

figure;
loglog(fpr, tpr);
xlabel('Probability of false alarm');
ylabel('Probability of detection');


% signalIndex = 1:50;
% noSignalIndex = 60:273;
% 
% x = mean(GoFBlockPd(signalIndex));
% y = mean(GoFBlockPd(noSignalIndex));

rmpath(['.' filesep 'data']);
