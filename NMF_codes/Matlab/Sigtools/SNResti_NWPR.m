function SNRestimate = SNResti_NWPR(Signal, M)

SignalBuff = buffer(Signal, M);

Pn = zeros(size(SignalBuff, 2), 1);
Pw = zeros(size(SignalBuff, 2), 1);

for i = 1:size(SignalBuff, 2)
    Pn(i) = sum(real(SignalBuff(:,i)))^2 + sum(imag(SignalBuff(:,i)))^2;
    Pw(i) = sum(real(SignalBuff(:,i)).^2 + imag(SignalBuff(:,i)).^2);
end
averagePn = mean(Pn);
averagePw = mean(Pw);

SNRestimate = pow2db((M*(averagePn/averagePw) - 1) / ((M - (averagePn/averagePw))));