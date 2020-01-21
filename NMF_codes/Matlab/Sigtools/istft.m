function [X,T] = istft(S,varargin)
%ISTFT Inverse short-time Fourier transform.
%   X = ISTFT(S) returns the inverse short-time Fourier transform (ISTFT)
%   of S. S must be a matrix where the rows correspond to frequency and the
%   columns correspond to time. S is expected to be two-sided and centered.
%
%   X = ISTFT(S,Fs) specifies the sample rate of X in hertz as a positive
%   scalar.
%
%   X = ISTFT(S,Ts) specifies Ts as a positive scalar duration
%   corresponding to the sample time of X. The sample rate in this case is
%   calculated as 1/Ts.
%
%   X = ISTFT(...,'Window',WINDOW) specifies the window used in calculating
%   the ISTFT. Perfect time-domain reconstruction requires the ISTFT window
%   to match the window used to generate the STFT. Use the function ISCOLA
%   to check a window/overlap combination for constant overlap-add (COLA)
%   compliance. COLA compliance is a requirement for perfect reconstruction
%   for non-modified spectra. The default is a Hann window of length 128.
%
%   X = ISTFT(...,'OverlapLength',NOVERLAP) specifies an integer number of
%   samples of overlap between adjoining segments. NOVERLAP must be smaller
%   than the length of the window. Perfect time-domain reconstruction
%   requires the ISTFT NOVERLAP to match the NOVERLAP used to generate the
%   STFT. Use the function ISCOLA to check a window/overlap combination for
%   constant overlap-add (COLA) compliance. COLA compliance is a
%   requirement for perfect reconstruction for non-modified spectra. If
%   NOVERLAP is not specified, it is set to the largest integer less than
%   or equal to 75% of the window length.
%
%   X = ISTFT(...,'FFTLength',NFFT) specifies the integer number of
%   frequency points used to calculate the discrete Fourier transform.
%   Perfect time-domain reconstruction requires the ISTFT NFFT to match the
%   NFFT used to generate the STFT. NFFT defaults to the length of WINDOW.
%
%   X = ISTFT(...,'Method',METHOD) specifies the inversion method to use:
%       'ola'  - Overlap-Add 
%       'wola' - Weighted Overlap-Add
%   If the method is set to 'wola', a second window is applied after the
%   inverse DFT and prior to the final overlap-add stage, and a Griffin-Lim
%   normalization is performed. Typically, the analysis window used to
%   generate the STFT is the same as the synthesis window used during the
%   ISTFT. METHOD defaults to 'wola'.
% 
%   X = ISTFT(...,'ConjugateSymmetric',CONJUGATESYMMETRIC) is specified as
%   a logical true if S is symmetric or a logical false if S is
%   nonsymmetric. When S (the STFT matrix input to the ISTFT function) is
%   not exactly conjugate symmetric due to round-off error,
%   CONJUGATESYMMETRIC set to true ensures S is treated as if it were
%   conjugate symmetric. If S is conjugate symmetric, then the inverse
%   transform computation is faster, and the output is real. This
%   name-value pair is not supported for code generation. Defaults to
%   false.
%
%   X = ISTFT(...,'Centered',CENTERED) treats the input S as a two-sided,
%   centered transform if CENTERED is true. The function rearranges S so
%   that the zero-frequency component is the first row of the array. If
%   CENTERED is false, S is not rearranged. CENTERED defaults to true.
%
%   [X,T] = ISTFT(...) returns time vector T. If a sample rate is provided,
%   T is a vector of time values in seconds. If a sample time is provided,
%   then T is a duration array with the same time format as the input. If
%   no time information is provided, the output is a vector of sample
%   numbers. X contains the reconstructed time-domain signal.
%
%    % EXAMPLE 1: 
%       % Compute the ISTFT of a real signal using the overlap-add method.
%       fs = 10240;
%       t = 0:1/fs:0.5-1/fs;
%       x = 5*sin(2*pi*t*10);
%       D = duration(0,0,1/fs);
%       win = hamming(512,'periodic');
%       S = stft(x,D,'Window',win,'OverlapLength',384,'FFTLength',1024);
%       [X,T] = istft(S,D,'Window',win,'OverlapLength',384,...
%           'FFTLength',1024,'Method','ola','ConjugateSymmetric',true); 
%       
%       % Plot original and resynthesized signals.  
%       plot(t,x,seconds(T),X,'-.')
%       axis tight
%       xlabel('Time (s)')
%       ylabel('Amplitude (V)')
%       title('Original and Reconstructed Signal')
%       legend('Original','Reconstructed')
% 
%    % EXAMPLE 2: 
%       % Compute the ISTFT of a complex signal using weighted overlap-add
%       % method.
%       fs = 3000;
%       t = 0:1/fs:1-1/fs;
%       x = exp(2j*pi*100*cos(2*pi*2*t))+randn(size(t))/100;
%       nwin = 100; 
%       win = hann(nwin,'periodic');
%       xZero = [zeros(1,nwin) x  zeros(1,nwin)]; % Zero pad to fix edges 
%       S = stft(xZero,fs,'Window',win,'OverlapLength',75);
%       [X,T] = istft(S,fs,'Window',win,'OverlapLength',75);
%       
%       % Remove zeros.
%       X(1:nwin) = []; 
%       X(end-nwin+1:end) = []; 
%       T = T(1:end-2*nwin); 
%       
%       % Plot original and resynthesized signals. 
%       plot(t,abs(x),T,abs(X),'-.')
%       xlabel('Time (s)')
%       ylabel('Amplitude (V)')
%       title('Original and Reconstructed Signal')
%       legend('Original','Reconstructed')    
%
%   See also STFT, ISCOLA, PSPECTRUM, IFFT

% [1] Crochiere, R. E. "A Weighted Overlap-Add Method of Short-Time
%     Fourier Analysis/Synthesis." IEEE Transactions on Acoustics, Speech,
%     and Signal Processing. Vol. ASSP-28, Feb. 1980, pp. 99-102.
% [2] Griffin, D. W. and J. S. Lim. "Signal Estimation from Modified
%     Short-Time Fourier Transform." IEEE Transactions on Acoustics,
%     Speech, and Signal Processing. Vol. ASSP-32, No. 2, April 1984.
% [3] Portnoff, M. R. "Time-Frequency Representation of Digital Signals
%     and Systems Based on Short-Time Fourier Analysis." IEEE Transactions
%     on Acoustics, Speech, and Signal Processing. Vol. ASSP-28, Feb. 1980,
%     pp. 55-69.
        
%   Copyright 2018-2019 The MathWorks, Inc.
%#codegen

%---------------------------------
% Check inputs/outputs
narginchk(1,14);
if coder.target('MATLAB') % For MATLAB
    nargoutchk(0,2);
else % For everything else 
    nargoutchk(1,2);
end

%---------------------------------
% Parse inputs
[data,opts] = signal.internal.stft.stftParser('istft',S,varargin{:});

%---------------------------------
% Compute ISTFT
[X,T] = computeISTFT(data,opts); 

%---------------------------------
% Update time vector if a duration was provided
if coder.target('MATLAB') && strcmpi(opts.TimeMode,'ts')
    T = seconds(T);
    T.Format = opts.TimeUnits;
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Helper functions
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [X,T] = computeISTFT(s,opts)
% Computes inverse short-time Fourier transform 
classCast = class(s); 

% Set variables 
win = opts.Window;
nwin = opts.WindowLength;
noverlap = opts.OverlapLength;
nfft = opts.FFTLength;
hop = nwin-noverlap;
nseg = opts.TimeAxisLength;
xlen = nwin + (nseg-1)*hop;
Fs = opts.EffectiveFs;

% Uncenter
if opts.Centered
    s = uncenter(s);
end

% IDFT
if coder.target('MATLAB') && opts.ConjugateSymmetric
    xifft = ifft(s,nfft,1,'symmetric');
else
    xifft = ifft(s,nfft,1,'nonsymmetric');
end

xifft = xifft(1:min(nwin,size(xifft,1)),:);

% Initialize time-domain signal
if isreal(xifft)
    x = zeros(xlen,1,classCast);
else
    x = complex(zeros(xlen,1,classCast));
end

% Set method
if strcmpi(opts.Method,'ola') 
    a = 0;
else % Else WOLA
    a = 1;
end

% Initialize normalization value
normVal = zeros(xlen,1);

% Overlap-add
for ii = 1:nseg
    x(((ii-1)*hop+1):((ii-1)*hop+nwin)) = x(((ii-1)*hop+1):((ii-1)*hop+nwin)) + xifft(:,ii).*win.^a;
    normVal(((ii-1)*hop+1):((ii-1)*hop+nwin)) = normVal(((ii-1)*hop+1):((ii-1)*hop+nwin))+win.^(a+1);
end

% Normalize
normVal(normVal<(nseg*eps)) = 1; % Don't normalize by small values
X = x./normVal;

% Time vector
T = cast(((0:numel(x)-1).')./Fs,classCast);

% Scale time vector in the case of normalized frequency
if opts.IsNormalizedFreq
    T = T.*opts.EffectiveFs; % sample
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function s = uncenter(s)
% Uncenter s

n = size(s,1);
if n/2==round(n/2)
  % even (nyquist is at end of spectrum)
  s = circshift(s,-(n/2-1));
else
  % odd
  s = ifftshift(s,1);
end
