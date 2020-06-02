function [iflawestiTF, IFrate, linearflag] = IFEstimator_pai(signal, iflaw, PerfIFflag)
% global Emuindex;
% global iflaw;               % real IF law
% global Intenumb;            % number of samples within integration time
% global HilbSignal;          % aggregated signal
% global lagwinlen;
% global Segnumb;             % segment number
% global Nonesegment;         % sampling number within one segment
% global PerfIFflag;
% global Krate;
% global iflawestiTF          % IF estimates based on TF analysis
% global phiesti;
% global IFrate;
% global iflawestiend;

Segnumb = 32;
W = 24;
Nonesegment = 1024;
Intenumb = length(signal);
Krate = 10; 
lagwinlen = 41;

% IF estimation
if PerfIFflag == 1
    iflawestiTF = iflaw;
else
    iflawestiTF = zeros(Intenumb,1);
    iflawestiTFform = zeros(Intenumb,1);
    for p = 1:Segnumb    
        % extract this segment signal
        HilbSignalpart = signal((p-1)*Nonesegment+1:p*Nonesegment);
        % perform the Modified-B distribution
        tfBD = quadtfd(HilbSignalpart,lagwinlen,1,'mb',0.05,Nonesegment);
        % IF estimation using MBD results
        [Y,I] = max(abs(tfBD),[],1);
        iflawestiTFform((p-1)*Nonesegment+1:p*Nonesegment) = (I-1)*0.5/Nonesegment; 
        iflawestiTF((p-1)*Nonesegment+1:p*Nonesegment) = I*0.5/Nonesegment; 
    end
end

%% IP estimation
% global phiestiend;
% phiesti = 2*pi*cumsum(iflawestiTF); 
% % phiesti = 2*pi*cumsum(iflawestiTF)+phiestiend; 
% phiesti = mod(phiesti,2*pi);
% 
% global phiestiform;
% phiestiform = 2*pi*cumsum(iflawestiTFform); 
% phiestiform = mod(phiestiform,2*pi);
% phiestiform = 0;
% phiestiend = phiesti(end);
%% Chirp rate estimation
IFrate = zeros(Intenumb,1);
for m = 1:Intenumb
    range = max(1,m-Krate):min(m+Krate,Intenumb);
    H = [(range)',iflawestiTF(range)];   
    C = cov(H);
    lamda = eig(C);
    IFrate(m) = (max(lamda)-C(1,1))/C(1,2); 
    
    if isnan(IFrate(m))
        Krate2 = 50;
        range = max(1,m-Krate2):min(m+Krate2,Intenumb);        
        H = [(range)',iflawestiTF(range)];  
        C = cov(H);
        lamda = eig(C);
        IFrate(m) = (max(lamda)-C(1,1))/C(1,2);        
    end
    
end

[IFratepeaks] = localpeaksforN(abs(IFrate),101,0.5); 
% for q = 1:Segnumb-1
%     if length(find(IFratepeaks==q*Nonesegment)) == 0
%         IFratepeaks = [IFratepeaks,q*Nonesegment];
%     end
% end

[IFratepeaks,IFratepeaksI] = sort(IFratepeaks,'ascend');
% global linearflag;
linearflag = ones(1,Intenumb);
% global W;
for i = 1:length(IFratepeaks)
    range = max(1,IFratepeaks(i)-W):min(Intenumb,IFratepeaks(i)+W-1);
    linearflag(range) = zeros(1,length(range));
end
   
linearflag(1:W) = zeros(1,W);
linearflag(Intenumb-W+1:Intenumb) = zeros(1,W);