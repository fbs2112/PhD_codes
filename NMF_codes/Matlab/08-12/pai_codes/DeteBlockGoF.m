%% Interference detection algorithm based on GoF test using block-wise STFT
function DeteBlockGoF
%% Global variables
global Emuindex;                                   % emulate index
global Loopnumb;                                   % loop number
global Segnumb;                                    % number of TF observation intervals within an integration time
global Nonesegment;                                % number of samples within a TF observation interval
global Signal;                                     % aggregated signal
global WinLBlock;                                  % window length used in block-wise STFT
global MBlock;                                     % number of samples for each block-wise STFT
global PfaVector
global GoFBlockDeteflag;                           % detection flag for GoF-based interference detection algorithm using block-wise STFT
global JNRIndex;
global periodIndex
global bandwidthIndex
global bandwithLoop
%% GoF-based interference detection algorithm using block-wise STFT
% Window function
h = window('rectwin',WinLBlock);

for index = 1:Segnumb
    % Compute block-wise STFT
    [tfr,~,~] = tfrstftblock(Signal((index-1)*Nonesegment+1:index*Nonesegment),1:MBlock,MBlock,h);
    % Apply GoF test to each frequency slice of block-wise STFT, and determine the detection flag
    for k = 1:MBlock
        [~, GoFBlockpvalue] = chi2gof((abs(tfr(k,:))).^2,'cdf',{@chi2cdf,2},'nparams',0);
        for pfaIndex = 1:length(PfaVector)
            
            if GoFBlockpvalue < PfaVector(pfaIndex)
                GoFBlockDeteflag(bandwidthIndex, periodIndex, JNRIndex, pfaIndex, k,(Emuindex-1)*Segnumb+index) = 1;
            else
                GoFBlockDeteflag(bandwidthIndex, periodIndex, JNRIndex, pfaIndex, k,(Emuindex-1)*Segnumb+index) = 0;
            end
        end
    end
end
% Compute the detection probability
if bandwidthIndex == bandwithLoop
    save(['.' filesep 'data' filesep 'resultsPai05.mat'], 'GoFBlockDeteflag');
end
