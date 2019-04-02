clear;
clc;
close all;

X = [1,2,3;4,5,6;7,8,9];
W0 = [1 1;8 2;1 3];
H0 = [1 2 2;3 3 8];


% opt = statset('Maxiter',1000,'Display','final');
% [W2, H2] = nmf_py(X, 2, 'w0', W, 'h0', H, 'options', opt, 'algorithm','mult');
% 
% numW = (X./(W*H + eps(X)))*H.';
% denW = (ones(size(H, 2))*H.') + eps(numW);
% 
% W = W .* (numW ./ denW);
% 
% numH = W.'*(X./(W*H + eps(X)));
% denH = W.'*ones(size(W, 1)) + eps(numH);
% 
% H = H .* (numH ./ denH);


X = [1,2,3;4,5,6;7,8,9];
W0 = [1 1 8].';
H0 = [1 2 2];

numW = (X./(W0*H0 + eps(X)))*H0.';
denW = (ones(size(W0, 1), size(H0, 2))*H0.') + eps(numW);

W = W0 .* (numW ./ denW); 

numH = W.'*(X./(W*H0 + eps(X)));
denH = W.'*ones(size(W, 1), size(H0, 2)) + eps(numH);

H = H0 .* (numH ./ denH);


%Do more tests with this implementation. It seems to be working



