clear;
clc;
close all;

frequency = '1152';
% fileName = '2018-09-01-09_52_55_0000020535312384.000000'; %train
fileName = '2018-09-01-09_52_55_0000000000000000.000000'; %test

load(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep fileName '_labels.mat']);
% save(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep fileName '_labels_bkp.mat'], 'interferenceDetFlag');

falseLabels = find(interferenceDetFlag == 0 | interferenceDetFlag == 1 | interferenceDetFlag == 4);
trueLabels = find(interferenceDetFlag == 2 | interferenceDetFlag == 3);

aux = [8 30 50 52 69 79 80 81 89 91 94 117 126 130 138 141 142 150 174 203 222 227 230 242 289 308 ...
    314 356 362 363 364 365 366 370 375 406 424 456 482 504 505 507 508 514 516 530 582 665 666 667 ...
    668 669 672 676 728 755 768 782 802 804 805 807 815 838 849 864 890];

interferenceDetFlag(falseLabels(aux)) = 3;
% save(['.' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep fileName '_labels.mat'], 'interferenceDetFlag');
