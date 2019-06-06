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
params.JNRVector = [10];

numberOfTrainingCells = 1000;
numberOfGuardCells = 100;
rng(random_state);

%signal mixture definition---------
t = 0:1/fs:(secondsOfData - 1/fs);
f = ((bandwidth/2)/secondsOfData)*t + f0;
f1 = 1;

signal1 = exp(1j*2*pi*f1*f.*t).';
% signal1 = repmat(signal1, ceil(100e-6/secondsOfData), 1);
signal1 = repmat(signal1, 2, 1);

onsetTime = 1000e-6;
offsetTime = onsetTime + length(signal1)/fs;
signal1 = [zeros(round(onsetTime*fs), 1); signal1; zeros(round(1000e-6*fs), 1)];
signal1Length = length(signal1);

onset = find(signal1, 1, 'first') - params.nperseg;
offset = signal1Length;

mixtureSignal = signal1;
%--------------------------------------------


%---------TPSW parameters-------------------
N = 300;
M = 5;
alpha = 4;
%-------------------------------------------

window_length = round(3e-6*fs);
window_median_length = 401;
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

detector = phased.CFARDetector('NumTrainingCells', numberOfTrainingCells, 'NumGuardCells', numberOfGuardCells,...
    'ProbabilityFalseAlarm', 0.05, 'ThresholdOutputPort', true);
detector.Method = 'SOCA';
% detector.Rank = 100;

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
            
            output = abs(output).^2;
            
            
%             outputVar = window_eval(output, window_length, @var);
%             outputVarTPSW = tpsw_filt(outputVar, M, N, alpha);
%             outputVar = outputVar ./ max(outputVar);
%             outputVarTPSW = outputVarTPSW ./ max(outputVarTPSW);
            
            %                 [detection_res, thres] = detector(outputVarTPSW.', 1:length(outputVarTPSW));
            %                 detection_res = window_eval(detection_res, window_median_length, @median);
            %                 figure;
            %                 plot(t*1e6, outputVar);
            %                 hold on;
            %                 plot(t*1e6, outputVarTPSW);
            %                 plot(t*1e6, detection_res);
            %                 plot(t*1e6, thres, '-.');
            %                 ylabel('Normalized Magnitude');
            %                 xlabel('Time [$\mu$s]');
            %                 ylim([0 1.1])
            %                 xlim([min(t) max(t)]*1e6);
            %
            %Drawing lines------------------------------------------------------
            %                 line([t(onset) t(onset)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
            %                 line([t(offset) t(offset)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
            
            %                 line([t(onset) t(onset)]*1e6, [0 1.1], 'Color','black','LineStyle','--');
            %                 line([t(onset  + round(5*window_length/2)) t(onset + round(5*window_length/2))]*1e6, [0 1.1], 'Color','black','LineStyle','--');
            %
            %                 h = findobj(gca,'Type','line');
            %                 legend([h(6) h(5) h(4)  h(3)], 'Variance', 'TPSW', 'Detector', 'Threshold')
            %
            %                 linewidth = 1.5;
            %                 fontname = 'Times';
            %                 fontsize = 24;
            %
            %                 figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
            %
            %                 dataPath = ['..' filesep 'figs' filesep '05-23' filesep];
            %
            %                 formatFig(gcf, [dataPath 'cfar_tpsw' num2str(params.JNRVector(JNRIndex))], 'en', figProp);
            
            %                 line([t(offset) t(offset)]*1e6, [0 1.1], 'Color','black','LineStyle','--');
            %                 line([t(offset + round(5*window_length/2)) t(offset + round(5*window_length/2))]*1e6, [0 1.1], 'Color','black','LineStyle','--');
            
            %                 line([t(onset) t(onset)]*1e6, [0 1.1], 'Color','magenta','LineStyle','--');
            %                 line([t(offset + round(window_length*2)) t(offset + round(window_length*2))]*1e6, [0 1.1], 'Color','magenta','LineStyle','--');
            %-------------------------------------------------------------------
            for thresholdIndex = 1:length(thresholdVector)
                %Detection assessment---------------------------
                [detection_res, threshold] = detector(output, 1:length(output));
                detection_res = window_eval(detection_res, window_median_length, @median);
                
                figure;
                plot(t*1e6, output)
                
                hold on;
                plot(t*1e6, detection_res);
                plot(t*1e6, threshold)
                legend('$\mathbf{s}$', ' $\mathbf{o}_{\mathrm{med}}$', 'CFAR Threshold');  
                xlim([0 max(t)*1e6]);
                xlabel('Time [$\mu$s]');
                ylabel('Magnitude');
                
                if any(detection_res(1:onset))
                    fp(loopIndex, JNRIndex, stdIndex,thresholdIndex) = fp(loopIndex, JNRIndex, stdIndex, thresholdIndex) + 1;
                else
                    tn(loopIndex, JNRIndex, stdIndex, thresholdIndex) = tn(loopIndex, JNRIndex, stdIndex, thresholdIndex) + 1;
                end
                
                if any(detection_res(onset+1:onset + round(5*window_length/2)))
                    tp(loopIndex, JNRIndex, stdIndex, thresholdIndex) = tp(loopIndex, JNRIndex, stdIndex, thresholdIndex) + 1;
                else
                    fn(loopIndex, JNRIndex, stdIndex, thresholdIndex) = fn(loopIndex, JNRIndex, stdIndex, thresholdIndex) + 1;
                end
                
                if any(detection_res(onset  + round(5*window_length/2) + 1:end))
                    fp(loopIndex, JNRIndex, stdIndex, thresholdIndex) = fp(loopIndex, JNRIndex, stdIndex, thresholdIndex) + 1;
                else
                    tn(loopIndex, JNRIndex, stdIndex, thresholdIndex) = tn(loopIndex, JNRIndex, stdIndex, thresholdIndex) + 1;
                end
                
                %-------------------------------------------------
            end
        end
    end
end

% save(['..' filesep '.' filesep 'data' filesep '05-23' filesep 'results11.mat'], 'fp', 'tp', 'tn', 'fn');

rmpath(['..' filesep '.' filesep 'Misc'])
rmpath(['..' filesep '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep '.' filesep 'Sigtools' filesep])