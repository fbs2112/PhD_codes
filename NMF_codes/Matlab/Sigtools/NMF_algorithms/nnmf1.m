function [w,h,dnorm] = nnmf1(a,w0,h0,ismult,maxiter,tolfun,tolx,...
                                   dispnum,repnum,usePool)
% Single non-negative matrix factorization
nm = numel(a);
sqrteps = sqrt(eps);


% Display progress.  For parallel computing, the replicate number will be
% displayed under the worker performing the replicate.
if dispnum>1 % 'final' or 'iter' 
    if usePool 
        labindx = internal.stats.parallel.workerGetValue('workerID');
        dispfmt = '%8d\t%8d\t%8d\t%14g\t%14g\n';
    else
        dispfmt = '%7d\t%8d\t%12g\t%12g\n';
    end   
end    
    
for j=1:maxiter
    if ismult
        % Multiplicative update formula
        numer = w0'*a;
        h = max(0,h0 .* (numer ./ ((w0'*w0)*h0 + eps(numer))));
        numer = a*h';
        w = max(0,w0 .* (numer ./ (w0*(h*h') + eps(numer))));
    else
        % Alternating least squares
        h = max(0, w0\a);
        w = max(0, a/h);
    end
    
    % Get norm of difference and max change in factors
    d = a - w*h;
    dnorm = sqrt(sum(sum(d.^2))/nm);
    dw = max(max(abs(w-w0) / (sqrteps+max(max(abs(w0))))));
    dh = max(max(abs(h-h0) / (sqrteps+max(max(abs(h0))))));
    delta = max(dw,dh);
    
    % Check for convergence
    if j>1
        if delta <= tolx
            break;
        elseif dnorm0-dnorm <= tolfun*max(1,dnorm0)
            break;
        elseif j==maxiter
            break
        end
    end

    if dispnum>2 % 'iter'
       if usePool 
           fprintf(dispfmt,labindx,repnum,j,dnorm,delta);
       else
           fprintf(dispfmt,repnum,j,dnorm,delta);
       end
    end

    % Remember previous iteration results
    dnorm0 = dnorm;
    w0 = w;
    h0 = h;
end

if dispnum>1   % 'final' or 'iter'
    if usePool 
       fprintf(dispfmt,labindx,repnum,j,dnorm,delta);
   else
       fprintf(dispfmt,repnum,j,dnorm,delta);
   end
end

end
