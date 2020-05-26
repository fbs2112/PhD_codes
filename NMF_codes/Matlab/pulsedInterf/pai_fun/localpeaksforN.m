%% Find local peaks of a sequence
function [peaks] = localpeaksforN(h,nhood,threcoef)

done = false;
hnew = h;
nhood_center = (nhood-1)/2;
peaks = [];

% Compute the threshold of the detected local peaks
h2 = h;
[maxY,maxI] = max(h2); 
thre = threcoef*max(h2);
% if maxI >= 70 && maxI <= length(h)-70   
%     thre = threcoef*max(h2);
% else
%     p1thre = maxI - nhood_center(1); p2thre = maxI + nhood_center(1);
%     h2(max(p1thre,1):min(p2thre,length(h2))) = 0;      
%     thre = threcoef*max(h2);
% end


while ~done
  [dummy,max_idx] = max(hnew); 
 
  if dummy > thre
    % Suppress this maximum and its close neighbors.
    p1 = max_idx - nhood_center(1); p2 = max_idx + nhood_center(1);
    hnew(max(p1,1):min(p2,length(h))) = 0;        
%     if max_idx >= 70 && max_idx <= length(h)-70       
    peaks = [peaks max_idx]; 
%     end
  else
    done = true;
  end
end
