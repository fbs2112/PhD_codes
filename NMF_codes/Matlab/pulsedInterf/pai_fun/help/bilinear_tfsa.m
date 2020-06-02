% GENERATE BILINEAR TIME-FREQUENCY TRANSFORMATIONS
%
% BILINEAR TRANSFORMATIONS transform the time domain signal (the
% variable INPUT SIGNAL) to an output quadratic time-frequency
% distribution (the variable OUTPUT TIME-FREQUENCY ARRAY)
%
% There are many different types of quadratic time-frequency distributions
% that can be specified:
%   
%   (i)    Wigner-Ville Distribution (WVD)
%   (ii)   Smoothed Wigner-Ville (SMOOTHED)
%   (iii)  Spectrogram (SPECX)
%   (iv)   Rihaczek-Margenau-Hill (RM)
%   (v)    Choi-Williams (CW)
%   (vi)   Born-Jordon (BJ)
%   (vii)  Zhao-Atles-Marks (ZAM)
%   (viii) Cross-Wigner-Ville Distribution (XWVD)
%   (ix)   B-distribution (B)
%   (x)    Modified B-distribution (MB)
%
%   and also an Ambiguity Function with transforms the signal into the
%   dopple-lag domain.
%
%   The various parameters associated with the distributions are as follows:
%
%   TIME-FREQUENCY ARRAY (tfd)
%
%      The computed time-frequency distribution.  size(tfd) will
%      return [a, b], where a is the next largest power of two above
%      FFT length, and b is floor(length(signal)/time_res) - 1.
%
%   INPUT SIGNAL
%
%      Input one dimensional signal to be analysed. An analytic signal
%      is required for this function, however, if signal is real, a
%      default analytic transformer routine will be called from this
%      function before computing tfd.
%
%   TIME RESOLUTION
%
%      The number of time samples to skip between successive slices.
%
%   LAG WINDOW LENGTH
%
%      This is the lag window length and controls the size of the
%      signal kernel (or instantaneous autocorrelation function) used
%      for analysis (lag_window_length must be odd). The kernel used
%      will be defined from -(lag_window_length+1)/2 to
%      +(lag_window_length+1)/2 in both time and lag dimensions.
%
%   FFT Length:
%
%      Zero-padding at the FFT stage of the analysis may be specified
%      by giving an FFT length larger than lag window length.  If FFT
%      length is not specified, or is smaller than the lag window
%      length, then the next highest power of two above lag window
%      length is used.  If FFT length is not a power of two, the next
%      highest power of two is used.
%
%    KERNEL OPTIONS
%
%      Each time-frequency distribution has various options that can
%      be adjusted to produce the desired tfd.
%
%
%
%  See Also:  quadtfd, analyt
