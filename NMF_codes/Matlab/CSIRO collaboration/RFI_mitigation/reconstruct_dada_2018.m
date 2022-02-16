clear;
clc;
close all;

frequency = '1152';
fileName = '2018-09-01-09_52_55_0000000000000000.000000';
resultsFileNumber = 'results6';
year = '2018';
filePath2 = ['E:' filesep 'Parkes data' filesep year filesep frequency filesep 'signal frames' filesep];

matObj = matfile(['..' filesep 'Labels' filesep year filesep frequency filesep fileName '_labels.mat']);
trueLabels = find(matObj.interferenceDetFlag == 1 | matObj.interferenceDetFlag == 3 | matObj.interferenceDetFlag == 5 | matObj.interferenceDetFlag == 6);
falseLabels = find(matObj.interferenceDetFlag == 0 | matObj.interferenceDetFlag == 2);

trueLabels = trueLabels(trueLabels <= 100);
falseLabels = falseLabels(falseLabels <= 100);

matResults = matfile(['.' filesep 'data' filesep year filesep frequency filesep resultsFileNumber '.mat']);
numberOfSignalFrames = length(falseLabels);

xA = squeeze(matResults.xHatA(:,2,1,:));
xB = squeeze(matResults.xHatB(:,2,1,:));

for i = 1:numberOfSignalFrames
    load([filePath2 fileName '_' num2str(falseLabels(i)) '.mat']);
    xA(:,falseLabels(i)) = parkesSignalA;
    xB(:,falseLabels(i)) = parkesSignalB;
end
clearvars parkesSignalA parkesSignalB;


% for i = 1:100
%     load([filePath2 fileName '_' num2str(i) '.mat']);
%     xA(:,i) = parkesSignalA;
%     xB(:,i) = parkesSignalB;
% end
% clearvars parkesSignalA parkesSignalB;

xA = xA(:);
xB = xB(:);

xA = buffer(xA, 2048, 0, 'nodelay');
xB = buffer(xB, 2048, 0, 'nodelay');

xBuffer = zeros(2048, size(xA, 2)*2);
xBuffer(:,1:2:end) = xA;
xBuffer(:,2:2:end) = xB;
xBuffer = xBuffer(:);

clearvars xA xB;
xBuffer2 = zeros(length(xBuffer)*2, 1);
xBuffer2(1:2:end,1) = real(xBuffer);
xBuffer2(2:2:end,1) = imag(xBuffer);
clearvars xBuffer;

xBuffer2 = xBuffer2 + 2^15;

filePath = ['E:' filesep 'Parkes data' filesep year filesep frequency filesep fileName '.dada'];
fileID = fopen(filePath, 'r', 'n', 'UTF-8');
header = fread(fileID, 4096, '*char');
fclose(fileID); 

filePath = ['E:' filesep 'Parkes data' filesep year filesep frequency filesep fileName '_' resultsFileNumber '.dada'];

% filePath = ['E:' filesep 'Parkes data' filesep year filesep frequency filesep fileName '_' 'testReconstruct' '.dada'];

fileID = fopen(filePath, 'w');
fwrite(fileID, header, '*char');    
fclose(fileID);

fileID = fopen(filePath, 'a');
fwrite(fileID, xBuffer2, 'uint16');
fclose(fileID);