function [code] = codeGen(numberOfCodes, codeLength)

% Generate PRN code
% Generate unipolar Gold code
Goldcode = zeros(numberOfCodes, codeLength);                   % unipolar code
G1 = ones(numberOfCodes,10);                            % register 1
G2 = ones(numberOfCodes, 10);                            % register 2
Goldphase = [2,6;3,7;4,8;1,8;2,9];                 % Gold code phase 

for index = 1:codeLength 
    for channelindex = 1:numberOfCodes
        Goldcode(channelindex,index) = mod(G1(channelindex,10)+mod(sum(G2(channelindex,Goldphase(channelindex,1))+G2(channelindex,Goldphase(channelindex,2))), 2), 2);
    end
    % feedback
    G1feedback = [G1(:,3) G1(:,10)];
    G2feedback = [G2(:,2) G2(:,3) G2(:,6) G2(:,8) G2(:,9) G2(:,10)];
    % shift
    G1(:,2:10) = G1(:,1:9);
    G2(:,2:10) = G2(:,1:9);
    G1(:,1) = mod(sum(G1feedback,2),2);
    G2(:,1) = mod(sum(G2feedback,2),2);
end

% Generate bipolar PRN code
code = 2*Goldcode-1;