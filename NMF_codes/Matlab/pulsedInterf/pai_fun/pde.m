% Finite phase difference and Kay frequency estimation
%
% Estimates the instantaneous frequency of the input signal using
% General Phase Difference (FFD, CFD, 4th, 6th) and Kay Frequency
% Estimation
%
% Usage:
%
%     ife = pde( signal, order [, window_length]);
%
% Inputs
%
%    signal
%
%	  Input one dimensional signal to be analysed.
%
%    order
%
%	  Order of the finite phase difference estimator. Available
%	  estimator orders are: 1, 2, 4, 6.
%
%    window_length
%
%	  Kay smoothing window length in the case of weighted phase
%	  difference estimator.  NB Set order = 1 for Kay's method.
%
%
%     For General Phase	Difference:
%
%     window_length:0;
%     order:1, 2, 4 or 6;
%
%     e.g. ife = pde( signal, order);
%
%     for Kay smoothing weighted phase difference:
%
%     window_length: desired window length;
%     order=1;
%
%     e.g. ife = pde( signal, 0, window_length );
%

% TFSA 5.5
% Copyright
% Signal Processing Research Laboratory
% Queensland University	of Technology
% GPO Box 2434
% Brisbane 4001
% Australia
%
% email: tfsa@qut.edu.au
