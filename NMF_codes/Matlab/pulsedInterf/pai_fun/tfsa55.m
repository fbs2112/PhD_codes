% tfsa5     Time-Frequency Signal Analysis toolbox
%
%
% The Signal Processing Research Centre grants you the right to
% install and use the enclosed software programs on a single
% computer. You may copy the software into any machine readable form
% for backup or archival purposes in support of your use of the
% Software on the single computer.  Plain text parts of the package
% may be customised for personal use, provided that any references to
% the authors and the copyright notices of the package are not
% removed. You may transfer the Software and license agreement to
% another party if the other party agrees to accept the terms and
% conditions of the Agreement, and if the Software remains unmodified
% by you.  If you transfer the Software, you must at the same time
% transfer all copies of the same and accompanying documentation, or
% destroy any copies not transfered. YOU WILL NOT: 1. Sublicense the
% Software; 2. Copy or transfer the Software in whole or in part,
% except as expressly provided for in the wording above;
% 3. Incorporate the Software in whole or in part into any commercial
% product.
%
% Although considerable effort has been expended to make the programs
% in TFSA 5.5 correct and reliable, we make no warranties, express or
% implied, that the programs contained in this package are free of
% error, or are consistent with any particular standard of
% merchantability, or that they will meet your requirements for any
% particular application.  They should not be relied on for solving a
% problem whose incorrect solution could result in injury to a person
% or loss of property.  If you do use the programs in such a manner,
% it is at your own risk.  The authors disclaim all liability for
% direct or consequential damages resulting from your use of this
% package.
%
% MATLAB is a trademark	of The MathWorks, Inc.
%
% UNIX is a trademark of American Telephone and Telegraph Company.
%
% MS-DOS is a trademark	of Microsoft Corporation.
%
% MS-WINDOWS is	a trademark of Microsoft Corporation.
%
% It is the intent of the Signal Processing Research Laboratory to
% continue to update TFSA 5.5 to reflect new ideas and algorithms,
% and to correct bugs which may be discovered.  Bug report are
% welcome, and should contain sufficient information to reliably
% reproduce the aberrant behaviour.  Telephone support is not
% available at present.  Suggested fixes or workarounds are welcome.
% Our address for matters concerning the TFSA 5.5 package is:
%
%     Email: tfsa@qut.edu.au
% or
%     Prof. B. Boashash
%     Signal Processing	Research Centre
%     School of	Electrical and Electronic Systems Engineering
%     Queensland University of Technology
%     GPO Box 2434
%     Brisbane 4001
%     Australia
%
% TFSA 5.5 was written by B. Boashash, B. Hatton, F. Zureiqat, P.
% Boles, A. Reilly, G. Frazer, B. Ristic, M. Roessgen, G. Roberts,
% J. Ralston., D. R. Iskander, B. Barkat and J. O' Toole.
%
% List of functions:
%
%   tfsa5       Graphical Interface to TFSA 5.5
%
%   analyt      Computes the analytic signal of	a real signal.
%   ambf        Computes the ambiguity function of input signal.
%   quadknl     Generate Time-Frequency Kernel Functions.
%   gsig        Generate various time- and frequency-varying test
%               signals.
%   lms         Adaptive least mean square instantaneous frequency
%               estimation.
%   pde         General order (1st, 2nd, 4th, 6th) and weighted phase
%               difference instantaneous frequency estimation.
%   pwvd4       Fourth order kernel polynomial Wigner-Ville distribution.
%   pwvd6       Sixth order kernel polynomial Wigner-Ville
%               distribution.
%   pwvpe       Peak of polynomial Wigner-Ville instantaneous frequency
%		estimation.
%   quadtfd     Generate quadratic Time-Frequency Representations.
%   rihazcek    Computes the Rihazcek Time-Frequency distribution
%               for a given signal.
%   rls         Adaptive recursive least squares instantaneous frequency
%               estimation.
%   sfpe        Peak of Spectrogram instantaneous frequency
%               estimation.
%   spec        Computes the Spectrogram Time-Frequency
%               distribution.
%   synthesize  Estimates time-domain signals from given
%               Time-Frequency distribution.
%   tfsapl      Time-frequency plot.
%   wfall       TFSA waterfall plot.
%   wlet        Foward and inverse Fast Wavelet Transforms.
%   wvd         Wigner-Ville distribution.
%   wvpe        Peak of Wigner-Ville instantaneous frequency estimation.
%   xwvd        Cross Wigner-Ville distribution.
%   zce         Zero-crossing instantaneous frequency estimation.
%
%
% Internal functions, not intended for direct use:
%
%   flatwf          Flat waterfall plot.
%   oploy           Function to find coefficients for polynomial.
%   tfsademo        DEMO frame.
%   tfsa_err_gui    Display error in frame.
%   tfsa_err        Function to handle errors displayed on command
%                   line.
%   tfsahelp        Help frame.
%   tfsamain        Main TFSA frame.
%   tfsamenu        Handles the menus in the main frame.
%   tfsaopen        Sets up main frame.
%   tfsa_plot2d     TFSA vector plot frame.
%   tfsa_wrn        Handles warnings displayed on command line.
%   uif_base        Template frame called by other uif_* frames.
%   uif_btfd        Callback to manage the bilinear TFD frame.
%   uif_defs        Declared constants used in uif_* frames.
%   uif_dirtfd      Direct method of implementation of some TFDs frame.
%   uif_gsig        Callback to manage the test signal generation frame.
%   uif_ife         Callback to manage the instantaneous frequency estimation
%                   frame.
%   uif_mtfd        Callback to manage the multi-linear tfd frame.
%   uif_plot        Callback to manage the plotting frame.
%   uif_synth       Callback to manage the synthesis frame.
%   uif_ts          Callback to manage the time-scale frame.
%   uideflts        Default values for all uicontrols.
%   unphase         Recovers phase of the analyic input signal.
%

% TFSA 5.5
% Copyright
% Signal Processing Research Centre
% Queensland University	of Technology
% GPO Box 2434
% Brisbane 4001
% Australia
%
% email: tfsa@qut.edu.au


function tfsa55()

tfsaopen('initialize');



