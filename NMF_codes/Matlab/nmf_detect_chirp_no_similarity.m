clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')


set(groot,'defaultFigureVisible','off')
% set(groot,'defaultFigureVisible','on')

addpath(['.' filesep 'Sigtools' filesep])
addpath(['.' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['.' filesep 'Misc'])

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;
figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);

dataPath = ['.' filesep 'figs' filesep '04-10' filesep];

fs = 32.768e6;
numberOfSources = 1;
secondsOfData = 8.62e-6;
secondsOfSilence = 100e-6;
numberOfSamples = secondsOfData*fs;
totalSamples = 4096;
bandwidth = 1e6;
f0 = 0;
random_state = 42;

params.fs = fs;
params.nfft = 128;
params.nperseg = 128;
params.overlap = params.nperseg-1;
params.hop_size = params.nperseg - params.overlap;
params.numberOfSources = numberOfSources;
params.init = 'random';
params.betaDivergence = 'kullback-leibler';
params.numberOfIterations = 10000;
params.tolChange = 1e-6;
params.tolError = 1e-6;
params.repetitions = 1;
params.JNRVector = [-10 -5 0];

rng(random_state);
save_fig = false;

%signal mixture definition---------
t = 0:1/fs:(secondsOfData - 1/fs);
f = ((bandwidth/2)/secondsOfData)*t + f0;
f1 = 1;

signal1 = exp(1j*2*pi*f1*f.*t).';
signal1 = repmat(signal1, ceil(100e-6/secondsOfData), 1);
onsetTime = 20e-6;
offsetTime = onsetTime + length(signal1)/fs;

signal1 = [zeros(round(onsetTime*fs), 1); signal1;zeros(round(onsetTime*fs), 1)];
signal1Length = length(signal1);

onset = find(signal1, 1, 'first') - params.nperseg;
offset = find(signal1, 1, 'last') - params.nperseg;

mixtureSignal = signal1;
%--------------------------------------------

window_length = round(3e-6*fs);
similarityName = 'inner';
expName = 'detect_chirp_no_similarity';
thresholdVector = 0.1:0.1:0.9;
% thresholdVector = 0.5;
monteCarloLoops = 50;

%Pre allocation------------------------------------------------------------
fp = zeros(monteCarloLoops, length(thresholdVector), length(params.JNRVector));
tp = zeros(monteCarloLoops, length(thresholdVector), length(params.JNRVector));
fn = zeros(monteCarloLoops, length(thresholdVector), length(params.JNRVector));
tn = zeros(monteCarloLoops, length(thresholdVector), length(params.JNRVector));
%--------------------------------------------------------------------------


for loopIndex = 1:monteCarloLoops
    loopIndex
    [W, H, reconstructError, PxxAux, f, t] = nmf_eval(mixtureSignal, params);
    
    for JNRIndex = 1:length(params.JNRVector)
        
        for thresholdIndex = 1:length(thresholdVector)
            
            HNorm = H{1, JNRIndex}(1,:) ./ max(H{1, JNRIndex}(1,:));
            figure;
            plot(t*1e6,  HNorm);
            ylabel('Normalized Magnitude');
            xlabel('Time [$\mu$s]');
            ylim([0 1.1])
            xlim([min(t) max(t)]*1e6);
            hold on
            
            %Drawing lines------------------------------------------------------
            line([t(onset) t(onset)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
            line([t(offset) t(offset)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
            
            line([t(onset  + round(window_length*2)) t(onset + round(window_length*2))]*1e6, [0 1.1], 'Color','black','LineStyle','--');
            line([t(offset - round(window_length*2)) t(offset - round(window_length*2))]*1e6, [0 1.1], 'Color','black','LineStyle','--');
            
            line([t(onset) t(onset)]*1e6, [0 1.1], 'Color','magenta','LineStyle','--');
            line([t(offset + round(window_length*2)) t(offset + round(window_length*2))]*1e6, [0 1.1], 'Color','magenta','LineStyle','--');
            %-------------------------------------------------------------------
            
            
            %Detection assessment---------------------------
            detection_res = detection_eval(HNorm, thresholdVector(thresholdIndex));
            plot(t.'*1e6, detection_res, '*-g')
            if any(detection_res(1:onset))
                fp(loopIndex, thresholdIndex, JNRIndex) = fp(loopIndex, thresholdIndex, JNRIndex) + 1;
            else
                tn(loopIndex, thresholdIndex, JNRIndex) = tn(loopIndex, thresholdIndex, JNRIndex) + 1;
            end
            
            if all(detection_res(onset  + round(window_length*2):offset - round(window_length*2)))
                tp(loopIndex, thresholdIndex, JNRIndex) = tp(loopIndex, thresholdIndex, JNRIndex) + 1;
            else
                fn(loopIndex, thresholdIndex, JNRIndex) = fn(loopIndex, thresholdIndex, JNRIndex) + 1;
            end
            
            if any(detection_res(offset + round(window_length*2):end))
                fp(loopIndex, thresholdIndex, JNRIndex) = fp(loopIndex, thresholdIndex, JNRIndex) + 1;
            else
                tn(loopIndex, thresholdIndex, JNRIndex) = tn(loopIndex, thresholdIndex, JNRIndex) + 1;
            end
            
            
            %------------------------------------------------
            
            if save_fig
                formatFig(gcf, [dataPath expName  '_' 'activation_'  num2str(params.JNRVector(JNRIndex))], 'en', figProp); %#ok<*UNRCH>
            end
            
        end
    end
end

save(['.' filesep 'data' filesep '04-10' filesep 'results04.mat'], 'fp', 'tp', 'tn', 'fn');

rmpath(['.' filesep 'Misc'])
rmpath(['.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['.' filesep 'Sigtools' filesep])