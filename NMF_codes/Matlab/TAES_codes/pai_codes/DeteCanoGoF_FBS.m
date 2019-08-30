function [pvalue, GoFBlockDeteflag] = DeteCanoGoF_FBS(signal, h, PfaVector)

GoFBlockDeteflag = zeros(length(signal), length(PfaVector));
pvalue = zeros(length(signal), 1);
[tfr,~,~] = tfrstft(signal, 1:length(signal), length(signal), h);

for k = 1:length(signal)
    [~, GoFBlockpvalue] = chi2gof((abs(tfr(k,:))).^2,'cdf',{@chi2cdf,2},'nparams',0);
    pvalue(k) = GoFBlockpvalue;
    for pfaIndex = 1:length(PfaVector)
        if GoFBlockpvalue < PfaVector(pfaIndex)
            GoFBlockDeteflag(k, pfaIndex) = 1;
        end
    end
end
