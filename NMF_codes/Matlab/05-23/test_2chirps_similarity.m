clear;
clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')

addpath(['..' filesep '.' filesep 'Sigtools' filesep])
addpath(['..' filesep '.' filesep 'Sigtools' filesep 'NMF_algorithms'])

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
% params.JNRVector = [-20 -15 -10 -5 0];
params.JNRVector = 0;


rng(random_state);

%signal mixture definition---------
t = 0:1/fs:(secondsOfData - 1/fs);
f = ((bandwidth/2)/secondsOfData)*t + f0;
f1 = 1;
f2 = 3;

signal1 = exp(1j*2*pi*f1*f.*t).';
signal1 = repmat(signal1, ceil(100e-6/secondsOfData), 1);
onsetTime1 = 20e-6;
offsetTime1 = onsetTime1 + length(signal1)/fs;

signal1 = [zeros(round(onsetTime1*fs), 1); signal1];
signal1Length = length(signal1);

onsetTime2 = 60e-6;

signal2 = exp(1j*2*pi*f2*f.*t).';
signal2 = repmat(signal2, ceil(100e-6/secondsOfData), 1);
signal2Length = length(signal2);
signal2 = [zeros(round(onsetTime2*fs), 1); signal2; zeros(round(onsetTime2*fs), 1)];


offsetTime2 = onsetTime2 + signal2Length/fs;

signal1 = [signal1;zeros(length(signal2) - signal1Length, 1)];

mixtureSignal = signal1 + signal2;

mixtureSignal = mixtureSignal(1:floor(offsetTime2*fs));

%--------------------------------------------

onset1 = find(signal1, 1, 'first') - params.nperseg;
offset1 = find(signal1, 1, 'last') - params.nperseg;

onset2 = find(signal2, 1, 'first') - params.nperseg;
offset2 = find(signal2, 1, 'last') - params.nperseg;


window_length = round(5e-6*fs);
window_median_length = 51; 
similarityName = 'inner';

stdVector = 0;

if strcmp(similarityName, 'gaussian')
    stdVector = 7:2:13;
end
expName = 'detect_chirp_similarity_inner_window_median';
thresholdVector = 0.1:0.1:0.9;
monteCarloLoops = 50;

%Pre allocation------------------------------------------------------------
fp_sim_window = zeros(monteCarloLoops, length(thresholdVector), length(params.JNRVector), length(stdVector));
tp_sim_window = zeros(monteCarloLoops, length(thresholdVector), length(params.JNRVector), length(stdVector));
fn_sim_window = zeros(monteCarloLoops, length(thresholdVector), length(params.JNRVector), length(stdVector));
tn_sim_window = zeros(monteCarloLoops, length(thresholdVector), length(params.JNRVector), length(stdVector));
%--------------------------------------------------------------------------


for loopIndex = 1:monteCarloLoops
    loopIndex
    [W, H, reconstructError, PxxAux, f, t] = nmf_eval(mixtureSignal, params);
    
    for JNRIndex = 1:length(params.JNRVector)
        
        for thresholdIndex = 1:length(thresholdVector)
            
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
                
                output(output<0) = 0;
                outputNorm = output ./ max(output);
                
                figure;
                plot(t*1e6, H{1, JNRIndex}(1,:)./ max(H{1, JNRIndex}(1,:)));
                hold on;
                line([t(onset1) t(onset1)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
                line([t(offset1) t(offset1)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
                line([t(onset2) t(onset2)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
                line([t(offset2) t(offset2)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
                title('H')
                
                figure;
                plot(t*1e6, outputNorm);
                hold on;
                line([t(onset1) t(onset1)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
                line([t(offset1) t(offset1)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
                line([t(onset2) t(onset2)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
                line([t(offset2) t(offset2)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
                title('sim')
                
                figure;
                hold on;
                outputMean = window_eval(output, window_length, @mean);
                plot(t*1e6, outputMean ./ max(outputMean));
                line([t(onset1) t(onset1)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
                line([t(offset1) t(offset1)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
                line([t(onset2) t(onset2)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
                line([t(offset2) t(offset2)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
                title('sim windowed')
                
                figure;
                hold on;
                outputVar = window_eval(H{1, JNRIndex}(1,:), window_length, @var);
                plot(t*1e6, outputVar ./ max(outputVar));
                line([t(onset1) t(onset1)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
                line([t(offset1) t(offset1)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
                line([t(onset2) t(onset2)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
                line([t(offset2) t(offset2)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
                title('H var')
                
                figure;
                hold on;
                outputVar = window_eval(output, window_length, @max_eval);
                plot(t*1e6, outputVar ./ max(outputVar));
                line([t(onset1) t(onset1)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
                line([t(offset1) t(offset1)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
                line([t(onset2) t(onset2)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
                line([t(offset2) t(offset2)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
                title('sim var')
                
                
                
%                 %Drawing lines------------------------------------------------------
%                 line([t(onset1) t(onset1)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
%                 line([t(offset1) t(offset1)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
%                 
%                 line([t(onset2) t(onset2)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
%                 line([t(offset2) t(offset2)]*1e6, [0 1.1], 'Color','red','LineStyle','--');
                
%                 line([t(onset  + round(window_length*2)) t(onset + round(window_length*2))]*1e6, [0 1.1], 'Color','black','LineStyle','--');
%                 line([t(offset - round(window_length*2)) t(offset - round(window_length*2))]*1e6, [0 1.1], 'Color','black','LineStyle','--');
%                 
%                 line([t(onset) t(onset)]*1e6, [0 1.1], 'Color','magenta','LineStyle','--');
%                 line([t(offset + round(window_length*2)) t(offset + round(window_length*2))]*1e6, [0 1.1], 'Color','magenta','LineStyle','--');
                %-------------------------------------------------------------------
                
                %Detection assessment---------------------------
                detection_res = detection_eval(outputNorm, thresholdVector(thresholdIndex), window_median_length);
                if any(detection_res(1:onset))
                    fp_sim_window(loopIndex, thresholdIndex, JNRIndex, stdIndex) = fp_sim_window(loopIndex, thresholdIndex, JNRIndex, stdIndex) + 1;
                else
                    tn_sim_window(loopIndex, thresholdIndex, JNRIndex, stdIndex) = tn_sim_window(loopIndex, thresholdIndex, JNRIndex, stdIndex) + 1;
                end
                
                if any(detection_res(onset  + round(window_length/2):onset  + round(3*window_length/2)))...
                        && any(detection_res(offset + round(window_length/2): offset + round(3*window_length/2)))
                    tp_sim_window(loopIndex, thresholdIndex, JNRIndex, stdIndex) = tp_sim_window(loopIndex, thresholdIndex, JNRIndex, stdIndex) + 1;
                else
                    fn_sim_window(loopIndex, thresholdIndex, JNRIndex, stdIndex) = fn_sim_window(loopIndex, thresholdIndex, JNRIndex, stdIndex) + 1;
                end
                
                if any(detection_res(onset  + round(3*window_length/2) + 1:offset + round(window_length/2)-1))
                    fp_sim_window(loopIndex, thresholdIndex, JNRIndex) = fp_sim_window(loopIndex, thresholdIndex, JNRIndex) + 1;
                else
                    tn_sim_window(loopIndex, thresholdIndex, JNRIndex) = tn_sim_window(loopIndex, thresholdIndex, JNRIndex) + 1;
                end
                
                if any(detection_res(offset + round(window_length*2):end))
                    fp_sim_window(loopIndex, thresholdIndex, JNRIndex, stdIndex) = fp_sim_window(loopIndex, thresholdIndex, JNRIndex, stdIndex) + 1;
                else
                    tn_sim_window(loopIndex, thresholdIndex, JNRIndex, stdIndex) = tn_sim_window(loopIndex, thresholdIndex, JNRIndex, stdIndex) + 1;
                end
                
                %-------------------------------------------------
            end
            
        end
    end
end

% save(['..' filesep '.' filesep 'data' filesep '05-02' filesep 'results06.mat'], 'fp_sim_window', 'tp_sim_window', 'tn_sim_window', 'fn_sim_window');

rmpath(['..' filesep '.' filesep 'Sigtools' filesep 'NMF_algorithms'])
rmpath(['..' filesep '.' filesep 'Sigtools' filesep])