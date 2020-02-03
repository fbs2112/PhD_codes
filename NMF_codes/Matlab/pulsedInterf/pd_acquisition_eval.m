clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');

addpath(['..' filesep 'signalsGeneration' filesep]);
addpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);
addpath(['..' filesep 'Misc'])
addpath(['..' filesep 'Sigtools' filesep])

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
dataPath = ['..' filesep 'figs' filesep '11-12' filesep];

load sim_params_2.mat;

load(['.' filesep 'data' filesep 'nmf_testing_13.mat']);

varN = 2;
Pfa = 0.05;
thresholdBorio = sqrt(-2*varN*log(Pfa));
fi = 0;
fs = paramsSignal.Freqsamp;
Nd = 81;
DopStep = 125; %(Hz)

trueDelay = 10e-6;
trueDoppler = 1.5e3;
paramsSignal.numberOfGPSSignals = 1;

[~, locC] = GPSGen(paramsSignal);   % Resample the code at data sampling frequency

JNRVector = [-15 -10 -5 10 30 50];
JNRVector = [0];

Tc = paramsSignal.Intetime;     % Coherent integration time in ms
Nc = floor(Tc * fs);

K = 1;  % Number of non-coherent integrations
codeDl = (0:(Nc - 1)) / fs; % Code delays

monteCarloLoops = 100;

% for nbitsIndex = 1:length(nbits)
for loopIndex = 1:monteCarloLoops
    loopIndex
%     locC = quantise_gps(locC, nbits(nbitsIndex));
    for JNRIndex = 1:length(JNRVector)
        
        xHatBorio = mixtureSignal(:,JNRIndex,loopIndex);
        zeroedSamples(JNRIndex,loopIndex) = length(find(abs(xHatBorio) > thresholdBorio));
        xHatBorio(abs(xHatBorio) > thresholdBorio) = 0;
        
        GPSSignals(1,:) = xHat(:,2,JNRIndex,loopIndex).';
        GPSSignals(2,:) = xHatBorio;
        GPSSignals(3,:) = mixtureSignal(:,JNRIndex,loopIndex).';
        
        for i = 1:3
            sspace = 0;                 % Search space were the results will be stored
            for ii = 1:K
                y = GPSSignals(i, (ii - 1) * Nc + (1:Nc) );   % use just 1 period of code at the time
                % Compute the search space for a single coherent integration epoch
                Tsspace = DftParallelCodePhaseAcquisition( y, locC, Nc, Nd, DopStep, fs, fi );
                sspace = sspace + Tsspace;  % Non-coherently accumulate the results
            end           
            
            [maxVal(i,JNRIndex,loopIndex), DopInd] = max(max(sspace.'));
            
        end
    end
end
averagePeakValue = mean(maxVal, 3);

rmpath(['..' filesep 'Sigtools' filesep])
rmpath(['..' filesep 'Misc'])
rmpath(['..' filesep 'signalsGeneration' filesep]);
rmpath(['..' filesep 'signalsGeneration' filesep 'sim_params']);