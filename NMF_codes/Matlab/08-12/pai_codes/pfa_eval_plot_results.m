clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')


addpath(['.' filesep 'data']);

load resultsPai07.mat;

PfaVector = logspace(-5, 0, 17);

for k = 1:2    
    for i = 1:length(PfaVector)
        x = squeeze(GoFBlockDeteflagCell{k, 1}(1,1,i,:,:));
        for j = 1:size(x, 2)
            fp(k,j,i) = sum(x(:,j));
            tn(k,j,i) = size(x(:,j), 1) - fp(k,j,i);
        end
    end

end
fpr = fp./(fp+tn);

averageFpr = squeeze(mean(fpr, 2));
stdFpr = squeeze(std(fpr, [], 2));

figure;
plot(PfaVector, averageFpr(1,:));
hold on;
plot(PfaVector, averageFpr(2,:));

xAux = linspace(0, 1, size(averageFpr, 2));
yAux = linspace(0, 1, size(averageFpr, 2));
plot(xAux, yAux, '--');
legend('$L = 3$', '$L = 19$', 'Random guess'); 
grid on;

figure;
loglog(PfaVector, averageFpr(1,:));
hold on;
loglog(PfaVector, averageFpr(2,:));
ylim([1e-5 1e0]);

xAux = logspace(-5, 0, size(averageFpr, 2));
yAux = logspace(-5, 0, size(averageFpr, 2));
plot(xAux, yAux, '--');
legend('$L = 3$', '$L = 19$', 'Random guess'); 
grid on;

rmpath(['.' filesep 'data']);
