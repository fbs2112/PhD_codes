function xHat = NF_pai(signal, iflawestiTF, kalpha, rho)

phiesti = 2*pi*cumsum(iflawestiTF);
phiesti = mod(phiesti,2*pi);

Intenumb = length(signal);
Inteesti = exp(-1i*phiesti);
Filterinput = Inteesti.*signal;
FilteroutputIIR = zeros(Intenumb,1);
lastoutput = 0;
lastinput = 0;
for i = 1:Intenumb
    FilteroutputIIR(i) = kalpha*rho*lastoutput + Filterinput(i) - rho*lastinput;
    lastoutput = FilteroutputIIR(i);
    lastinput = Filterinput(i);
end
xHat = (FilteroutputIIR.*conj(Inteesti));
