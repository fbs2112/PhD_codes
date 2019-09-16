function [GoFBlockDeteflag, pvalue, kvalue] = deteTimeSlice(signal, params, thresholdVector)

timeBins = floor((length(signal) - params.overlap)/(params.nperseg - params.overlap));
signalBuffered = buffer(signal, params.nperseg, params.overlap, 'nodelay');
signalBuffered = signalBuffered(:,1:timeBins);
signalBuffered = [signalBuffered;zeros(params.nfft - params.nperseg, timeBins)]; 
PxxAux = dftmtx(params.nfft)*signalBuffered;

% [PxxAux, ~, ~] = spectrogram(signal, params.window, params.overlap, params.nfft, params.fs, 'centered', params.specType);

GoFBlockDeteflag = zeros(size(PxxAux, 2), length(thresholdVector));
pvalue = zeros(size(PxxAux, 2), 1);
kvalue = zeros(size(PxxAux, 2), 1);

for k = 1:size(PxxAux, 2)
%     [~, pvalue(k)] = chi2gof(abs(PxxAux(:,k)).^2,'cdf',{@chi2cdf,2},'nparams',0);
    [~, pvalue(k), ~, ~] = lillietest(abs(PxxAux(:,k)).^2, 'Distribution', 'exponential', 'Alpha', 0.05, 'MCTol',1e-2);
    for thresholdIndex = 1:length(thresholdVector)
        if pvalue(k) < thresholdVector(thresholdIndex)
            GoFBlockDeteflag(k, thresholdIndex) = 1;
        end
    end
end