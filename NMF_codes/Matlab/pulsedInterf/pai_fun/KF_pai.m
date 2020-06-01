function xHat = KF_pai(signal, IFEstimationFlag, iflawestiTF, IFrate, linearflag)

Intenumb = length(signal);
intesti = zeros(Intenumb, 1);

H = [1,0];
Pesti = diag([100,100]);
Xesti = zeros(2,1);  

for m = 1:Intenumb
    F = [2*cos(2*pi*iflawestiTF(m)),-exp(-1i*2*pi*IFrate(m));1,0];
    
    if IFEstimationFlag == 1
        if linearflag(m) == 1
            Q = diag([0.001,0]);
            R = 1;
        elseif linearflag(m) == 0
            Q = diag([1,0]);
            R = 0.01;
        end
    else
        if linearflag(m) == 1
            Q = diag([0.02,0]);
            R = 1;
        elseif linearflag(m) == 0
            Q = diag([1,0]);
            R = 0.01;
        end
    end
    
    Z = signal(m);
    Xpred = F*Xesti;
    Ppred = F*Pesti*F'+Q;
    KFG = Ppred*H'/(H*Ppred*H'+R);
    Pesti = (eye(2)-KFG*H)*Ppred*(eye(2)-KFG*H)'+KFG*R*KFG';
    Xesti = Xpred+KFG*(Z-H*Xpred);
    intesti(m) = Xesti(1);
end

xHat = (signal-intesti);