function GoFBlockDeteflag = deteTimeSlice(signal, h, MBlock, PfaVector)

GoFBlockDeteflag = zeros(MBlock, length(PfaVector));
[tfr,~,~] = tfrstftblock(signal, 1:MBlock, 256, h);

for k = 1:MBlock
    [~, GoFBlockpvalue] = chi2gof(abs(tfr(:,k)).^2,'cdf',{@chi2cdf,2},'nparams',0);
    for pfaIndex = 1:length(PfaVector)
        if GoFBlockpvalue < PfaVector(pfaIndex)
            GoFBlockDeteflag(k, pfaIndex) = 1;
        end
    end
end

