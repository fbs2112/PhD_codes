function SNRestimate = SNResti_pai(Signal)
 
%     SNRestimateRetr(index) = (mean(real(PresultRetr)).^2+mean(imag(PresultRetr)).^2)./var(PresultRetr);
PowerRetr = (real(Signal)).^2 + (imag(Signal)).^2;
meanPower = mean(PowerRetr);
varPower = var(PowerRetr);
SNRestimate = sqrt(meanPower^2-varPower)/(meanPower-sqrt(meanPower^2-varPower));
SNRestimate = 10*log10(SNRestimate);
