clear;
clc;
close all;

Fs = 128e6;
frequency = '1152';
fileName = '2018-09-01-09_52_55_0000000000000000.000000';
diskLetter = 'E:';
year = '2018';

filePath = ['E:' filesep 'Parkes data' filesep year filesep frequency filesep fileName '.dada'];
filePath2 = ['E:' filesep 'Parkes data' filesep year filesep frequency filesep fileName '_testReconstruct.dada'];

dataLength = round(10e-3*Fs) * 4; %4 is the number of bits

fileID = fopen(filePath, 'r','n','UTF-8');
dat = fread(fileID, 4096, '*char');
% fseek(fileID, 4096, 'bof');

fileID2 = fopen(filePath2, 'r','n','UTF-8');
dat2 = fread(fileID2, 4096, '*char');
% fseek(fileID2, 4096, 'bof');

fprintf(dat2);
fprintf('\n');

counter = 1;
i = 1;
while counter
    dat = fread(fileID, [dataLength, 1], 'uint16');
    datAux = dat;
    if any(isempty(dat))
        break;
    end

    dat = dat - 2^15; % offset binary
    dat_c = dat(1:2:end) + 1j*dat(2:2:end);
    dat_c = buffer(dat_c, 2048, 0, 'nodelay');
    
    parkesSignalA = dat_c(:,1:2:end);
    parkesSignalB = dat_c(:,2:2:end);
    parkesSignalA = parkesSignalA(:);
    parkesSignalB = parkesSignalB(:);
    
    dat2 = fread(fileID2, [dataLength, 1], 'uint16');
    errorDat = isequal(datAux, dat2);
    if any(isempty(dat2))
        break;
    end

    dat2 = dat2 - 2^15; % offset binary
    dat_c2 = dat2(1:2:end) + 1j*dat2(2:2:end);
    dat_c2 = buffer(dat_c2, 2048, 0, 'nodelay');
    
    parkesSignalA2 = dat_c2(:,1:2:end);
    parkesSignalB2 = dat_c2(:,2:2:end);
    parkesSignalA2 = parkesSignalA2(:);
    parkesSignalB2 = parkesSignalB2(:);
        
    errorA = isequal(parkesSignalA, parkesSignalA2);  
    errorB = isequal(parkesSignalB, parkesSignalB2);  
    
%     if ~errorA || ~errorB
%         error('DADA file reconstruction error');
%     end
   
    i = i + 1;
end
fclose(fileID);
fclose(fileID2);