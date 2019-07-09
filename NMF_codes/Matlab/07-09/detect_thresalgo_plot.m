clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
addpath(['..' filesep '.' filesep 'Misc'])

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
params.JNRVector = [-20 -15 -10 -5 0];
params.JNRVector = [-10];


rng(random_state);

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
window_median_length = 51;
similarityName = 'inner';
stdVector = 0;

if strcmp(similarityName, 'gaussian')
    stdVector = 7:2:13;
end

thresholdVector = 1;
monteCarloLoops = 50;

%Pre allocation------------------------------------------------------------
fp = zeros(monteCarloLoops, length(params.JNRVector), length(stdVector), length(thresholdVector));
tp = zeros(monteCarloLoops, length(params.JNRVector), length(stdVector), length(thresholdVector));
fn = zeros(monteCarloLoops, length(params.JNRVector), length(stdVector), length(thresholdVector));
tn = zeros(monteCarloLoops, length(params.JNRVector), length(stdVector), length(thresholdVector));
%--------------------------------------------------------------------------

for loopIndex = 1:monteCarloLoops
    loopIndex
    [W, H, reconstructError, PxxAux, f, t] = nmf_eval(mixtureSignal, params);
    
    for JNRIndex = 1:length(params.JNRVector)
        
        inputNMF = abs(PxxAux{1, JNRIndex}).^2;
        output = zeros(length(t), 1);
        
        for stdIndex = 1:length(stdVector)
            
            for tIndex = 1:length(t)
                if strcmp(similarityName, 'inner')
                    output(tIndex) = similarity_eval(inputNMF(:,tIndex), W{1, JNRIndex}(:,1), similarityName, 'normalized', true);
                else
                    output(tIndex) = similarity_eval(inputNMF(:,tIndex), W{1, JNRIndex}(:,1), similarityName, stdVector(stdIndex), true);
                end
            end
            
            output = TK_filtering(output);            
            for thresholdIndex = 1:length(thresholdVector)
                %Detection assessment---------------------------
%                 [detection_res, threshold] = detector(output, 1:length(output));
                
                [detection_res, ~, ~] = ThresholdingAlgo(output, 400, 3, 0.5);
                x = output./max(output);
                figure;
                plot(t*1e6, x)
                hold on;
                plot(t*1e6, detection_res, 'r--');
%                 figure
%                 [a,b] = findpeaks(detection_res);
%                 plot(findpeaks(detection_res));
                
                line([t(onset) t(onset)]*1e6, [-0.3 1.1], 'Color','black','LineStyle','-.');
                line([t(offset + 100) t(offset + 100)]*1e6, [-0.3 1.1], 'Color','black','LineStyle','-.');
                
                 ylabel('Normalized Magnitude');
                 xlabel('Time [$\mu$s]');
                 ylim([-0.3 1.1])
                 xlim([min(t) max(t)]*1e6);
                 
                 linewidth = 1.5;
                 fontname = 'Times';
                 fontsize = 24;
                 
                 figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
                 
                 dataPath = ['..' filesep 'figs' filesep '07-09' filesep];
                 formatFig(gcf, [dataPath 'tk' num2str(params.JNRVector(JNRIndex))], 'en', figProp);
%                 figure
%                 z = x(1:end - window_length+1) - y(window_length:end).';
%                 plot(z);
%                 
                
                
                figure;
                plot(t*1e6, detection_res);
                ylim([0 1.1]);
                figure;
                plot(t*1e6, threshold)
%                 legend('$\mathbf{s}$', ' $\mathbf{o}_{\mathrm{med}}$', 'CFAR Threshold');  
                xlim([0 max(t)*1e6]);
%                 xlabel('Time [$\mu$s]');
%                 ylabel('Magnitude');
                

                figure;
                histogram(output, 50);
                
%                if any(detection_res(1:onset))
%                    fp(loopIndex, thresholdIndex, JNRIndex, stdIndex) = fp(loopIndex, thresholdIndex, JNRIndex, stdIndex) + 1;
%                else
%                    tn(loopIndex, thresholdIndex, JNRIndex, stdIndex) = tn(loopIndex, thresholdIndex, JNRIndex, stdIndex) + 1;
%                end
%                
%                if any(detection_res(onset+1:onset + round(5*window_length/2)))
%                    tp(loopIndex, thresholdIndex, JNRIndex, stdIndex) = tp(loopIndex, thresholdIndex, JNRIndex, stdIndex) + 1;
%                else
%                    fn(loopIndex, thresholdIndex, JNRIndex, stdIndex) = fn(loopIndex, thresholdIndex, JNRIndex, stdIndex) + 1;
%                end
%                
%                if any(detection_res(onset  + round(5*window_length/2) + 1:end))
%                    fp(loopIndex, thresholdIndex, JNRIndex) = fp(loopIndex, thresholdIndex, JNRIndex) + 1;
%                else
%                    tn(loopIndex, thresholdIndex, JNRIndex) = tn(loopIndex, thresholdIndex, JNRIndex) + 1;
%                end
%                
%                if any(detection_res(onset+1:onset + round(5*window_length/2)))
%                    tp(loopIndex, thresholdIndex, JNRIndex, stdIndex) = tp(loopIndex, thresholdIndex, JNRIndex, stdIndex) + 1;
%                else
%                    fn(loopIndex, thresholdIndex, JNRIndex, stdIndex) = fn(loopIndex, thresholdIndex, JNRIndex, stdIndex) + 1;
%                end
%                 
%                 %-------------------------------------------------
            end
        end
    end
end

% save(['..' filesep '.' filesep 'data' filesep '06-11' filesep 'results11.mat'], 'fp', 'tp', 'tn', 'fn');

rmpath(['..' filesep '.' filesep 'Misc'])
rmpath(['..' filesep '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep '.' filesep 'Sigtools' filesep])