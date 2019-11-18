function sspace = DftParallelCodePhaseAcquisition(sig, locC, N, Nd, DopStep, fs, fi)
% Purpose.
% Evaluate the search space for the code acquisition using the time domain FFT technique   
%
% Syntax.
% sspace = DftParallelCodePhaseAcquisition(sig, locC, N, dstep, DopStep)  
%
% Input Parameters.
% sig       : [vector] the Galileo/GPS input signal, corrupted by Doppler
%              shift, code delay, noise and eventually interferer
% locC      : [vector] the local code replica
%
% N         : [integer] the input signal and local code length
%
% Nd        : [integer] number of Doppler bin used for the search space
%
% DopStep   : [Hz] Doppler bin width in Hz
%
% fs        : [Hz] sampling frequency
%
% fi        : [Hz] intermediate frequency 
%
% Output Parameters.
% sspace       : [matrix] the search space
%
% Author.
% Daniele Borio
%
% Version.
%   6 - 3 - 2006

fif = fi/fs;    % normalized intermediate frequency
deltaf = DopStep/fs;    % normalized Doppler step

sspace = zeros(Nd, N);

F_CA = conj(fft(locC));% /N;
t = 0:(N - 1);  % time index

for ff = 1:Nd,  
    fc = fif + (ff - ceil(Nd/2))*deltaf;        
    IQ_comp = exp(-2*j*pi*fc.*t).*sig;
    X = fft(IQ_comp);
    sspace(ff,:) = ifft(X.*F_CA);
end

sspace = real(sspace.*conj(sspace));

