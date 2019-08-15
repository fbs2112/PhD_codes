%% Interference detection algorithm based on GoF test using canonical STFT
function DeteCanoGoF
%% Global variables
global Emuindex;                                   % emulate index
global Loopnumb;                                   % loop number
global Segnumb;                                    % number of TF observation intervals within an integration time
global Nonesegment;                                % number of samples within a TF observation interval
global WinLCano;                                   % window length used in canonical STFT
global WintypeCano;                                % window type used in canonical STFT
global Signal;                                     % aggregated signal
global PfaVector

global GoFCanoDeteflag;                            % detection flag for GoF-based interference detection algorithm using canonical STFT
%% GoF-based interference detection algorithm using canonical STFT
% Window function
if WintypeCano == 1
    h = window('rectwin',WinLCano);
elseif WintypeCano == 2
    h = window('hamming',WinLCano);
elseif WintypeCano == 3
    h = window('gausswin',WinLCano);
end

for index = 1:Segnumb
    % Compute canonical STFT
    [tfr,~,~] = tfrstft(Signal((index-1)*Nonesegment+1:index*Nonesegment),1:Nonesegment,Nonesegment,h);
    % Apply GoF test to each frequency slice of canonical STFT, and determine the detection flag
    for k = 1:Nonesegment
        [~,GoFCanopvalue] = chi2gof((abs(tfr(k,:))).^2,'cdf',{@chi2cdf,2},'nparams',0);
        for pfaIndex = 1:length(PfaVector)
            
            if GoFCanopvalue < PfaVector(pfaIndex)
                GoFCanoDeteflag(pfaIndex, k,(Emuindex-1)*Segnumb+index) = 1;
            else
                GoFCanoDeteflag(pfaIndex, k,(Emuindex-1)*Segnumb+index) = 0;
            end
        end
    end
    
end
% Compute the detection probability
if Emuindex == Loopnumb
    %         GoFBlockdete(:,pfaIndex) = sum(GoFBlockDeteflag,2)/(Loopnumb*Segnumb);
%     save(['.' filesep 'data' filesep 'resultsPai04_' num2str(JNRIndex) '.mat'], 'GoFCanoDeteflag');
    save(['.' filesep 'data' filesep 'resultsPai04_2.mat'], 'GoFCanoDeteflag');

end