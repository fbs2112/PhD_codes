clear;
clc;
close all;

warning off;

frequency = '1152';
fileName = '2018-09-01-09_52_55_0000000000000000.000000'; %test

load(['..' filesep 'CSIRO collaboration' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep fileName '_labels_train_test.mat']);

params.fs = 128e6;
params.nfft = 256;
params.nperseg = 256;
params.overlap = 0;
params.hop_size = params.nperseg - params.overlap;
params.window = ones(params.nperseg, 1);
params.JNRVector = 0;

% thresholdVector = 200:0.5:3000;

% window_median_length_vector = 0;

monteCarloLoops = length(trueLabelsTest);

WLin = 32;     % Internal word length (for Goertzel's Algorithm)
FLin = 15;     % Internal fractional length (for Goertzel's Algorithm)
WL = 16;       % Word length (input data)
FL = 0;        % Fractional length (input data)

% Define quantiser objects
q = quantizer('fixed', 'convergent', 'Saturate', [WL FL]);
q0 = quantizer('fixed', 'convergent', 'Saturate', [WLin FLin]);

% fiObj = fimath('OverflowAction', 'Saturate', 'RoundingMethod','convergent');

numTrainingCells = 32;
numGuardCells = 16;
cfar = phased.CFARDetector('NumTrainingCells', numTrainingCells,'NumGuardCells', numGuardCells);
cfar.ThresholdFactor = 'Auto';
cfar.ProbabilityFalseAlarm = 1e-3;

k = 133;

% detection_res = zeros(monteCarloLoops, length(thresholdVector), length(window_median_length_vector));
% feature = zeros(monteCarloLoops, 1);

output_1090 = zeros(monteCarloLoops, 5000);
% output_1089 = zeros(monteCarloLoops, 5000);
% output_1091 = zeros(monteCarloLoops, 5000);

coef = quantize(q0,2*cos((2*pi*(k-1))/params.nfft));
out_coef = quantize(q0,exp((-1j*(2*pi*(k-1))/params.nfft)));

for loopIndex = 1:monteCarloLoops
    loopIndex
    
    load(['..' filesep 'CSIRO collaboration' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep fileName '_' num2str(trueLabelsTest(loopIndex)) '.mat']);
    parkesSignalBufferA = buffer(parkesSignalA, params.nfft, params.overlap, 'nodelay');
    N = length(parkesSignalBufferA);
%     detAux = max(abs(goertzel(parkesSignalBufferA, k, 1)));
    parkesSignalBufferA_Q = quantize(q, parkesSignalBufferA);
%     parkesSignalBufferA_Q2 = fi(parkesSignalBufferA, 1, WL, FL, fiObj);
    
    %     % 1090.0MHz detector
    %     coef_unqu = 2*cos((2*pi*(k-1))/params.nfft);
    %     out_coef_unqu = exp((-1j*(2*pi*(k-1))/params.nfft));
    %     dft_data = zeros(N, 1);
    %     for l = 1:N
    %         v_unqu = zeros(params.nfft, 1);
    %         for m = 1:params.nfft
    %
    %             if m == 1
    %                 v_unqu(m) = parkesSignalBufferA(m,l);
    %             elseif m==2
    %                 v_unqu(m) = coef_unqu*v_unqu(m-1)+parkesSignalBufferA(m,l);
    %             else
    %                 v_unqu(m) = coef_unqu*v_unqu(m-1)-v_unqu(m-2)+parkesSignalBufferA(m,l);
    %             end
    %             if m==params.nfft
    %                 dft_data(l,1) = abs(-(v_unqu(params.nfft)-(out_coef_unqu*v_unqu(params.nfft-1))));
    %             end
    %         end
    %     end
    
    % 1090.0MHz quantised detector
    
    dft_data_1090_quant = zeros(N, 1);
    for l = 1:N
        
        v = zeros(params.nfft, 1);
        for m = 1:params.nfft
            
            if m==1
                v(m) = parkesSignalBufferA_Q(m,l);
            elseif m==2
                v(m) = quantize(q0,coef*v(m-1)+parkesSignalBufferA_Q(m,l));
            else
                v(m) = quantize(q0,coef*v(m-1)-v(m-2)+parkesSignalBufferA_Q(m,l));
            end
            if m==params.nfft
                dft_data_1090_quant(l,loopIndex) = quantize(q0,abs(-(v(params.nfft)-(out_coef*v(params.nfft-1)))).^2);
            end
        end
    end
    
%     % 1090.0MHz quantised detector
%     coef = fi(2*cos((2*pi*(k-1))/params.nfft), 1, WLin , FLin, fiObj);
% %     coef = quantize(q0,2*cos((2*pi*(k-1))/params.nfft));
%     out_coef = fi(exp((-1j*(2*pi*(k-1))/params.nfft)), 1, WLin, FLin, fiObj);
% %     out_coef = quantize(q0,exp((-1j*(2*pi*(k-1))/params.nfft)));
%     dft_data_1090_quant2 = zeros(N, 1);
%     for l = 1:N
%         l
%         v = zeros(params.nfft, 1);
%         for m = 1:params.nfft
%             
%             if m==1
% %                 v(m) = quantize(q,parkesSignalBufferA_Q(m,l));
%                 v(m) = parkesSignalBufferA_Q2(m,l);
%             elseif m==2
% %                 v(m) = quantize(q0,coef*v(m-1)+parkesSignalBufferA_Q(m,l));
%                 v(m) = fi(coef*v(m-1)+parkesSignalBufferA_Q2(m,l), 1, WLin , FLin, fiObj);
%             else
%                 v(m) = fi(coef*v(m-1)-v(m-2)+parkesSignalBufferA_Q2(m,l), 1, WLin , FLin, fiObj);
% %                 v(m) = quantize(q0,coef*v(m-1)-v(m-2)+parkesSignalBufferA_Q(m,l));
%             end
%             if m==params.nfft
%                 dft_data_1090_quant2(l,loopIndex) = fi(abs(-(v(params.nfft)-(out_coef*v(params.nfft-1)))), 1, WLin , FLin, fiObj);
% %                 dft_data_1090_quant2(l,loopIndex) = quantize(q0,abs(-(v(params.nfft)-(out_coef*v(params.nfft-1)))));
%             end
%         end
%     end
    
    %     % 1089.5 MHz quantised detector
    %     coef = quantize(q0,2*cos((2*pi*(k-2))/params.nfft));
    %     out_coef = quantize(q0,exp((-1j*(2*pi*(k-2))/params.nfft)));
    %     dft_data_1089_quant = zeros(N, 1);
    %     for l = 1:N
    %         v = zeros(params.nfft, 1);
    %        for m = 1:params.nfft
    %            if m==1
    %                v(m) = quantize(q,parkesSignalBufferA(m,l));
    %            elseif m==2
    %                v(m) = quantize(q0,coef*v(m-1)+parkesSignalBufferA(m,l));
    %            else
    %                v(m) = quantize(q0,coef*v(m-1)-v(m-2)+parkesSignalBufferA(m,l));
    %            end
    %            if m==params.nfft
    %                dft_data_1089_quant(l,loopIndex) = quantize(q0,abs(-(v(params.nfft)-(out_coef*v(params.nfft-1)))));
    %            end
    %        end
    %     end
    %
    %     % 1090.5 MHz quantised detector
    %     coef = quantize(q0,2*cos((2*pi*(k))/params.nfft));
    %     out_coef = quantize(q0,exp((-1j*(2*pi*(k))/params.nfft)));
    %     dft_data_1091_quant = zeros(N, 1);
    %     for l = 1:N
    %         v = zeros(params.nfft, 1);
    %        for m = 1:params.nfft
    %            if m==1
    %                v(m) = quantize(q,parkesSignalBufferA(m,l));
    %            elseif m==2
    %                v(m) = quantize(q0,coef*v(m-1)+parkesSignalBufferA(m,l));
    %            else
    %                v(m) = quantize(q0,coef*v(m-1)-v(m-2)+parkesSignalBufferA(m,l));
    %            end
    %            if m==params.nfft
    %                dft_data_1091_quant(l,loopIndex) = quantize(q0,abs(-(v(params.nfft)-(out_coef*v(params.nfft-1)))));
    %            end
    %        end
    %     end
    %
    %
    %
        output_1090(loopIndex,:) = cfar(dft_data_1090_quant(:,loopIndex), 1:length(dft_data_1090_quant(:,loopIndex)));
    %     output_1089(loopIndex,:) = cfar(dft_data,1:length(dft_data_1089_quant));
    %     output_1091(loopIndex,:) = cfar(dft_data,1:length(dft_data_1091_quant));
    
end

save(['..' filesep 'CSIRO collaboration' filesep 'data' filesep 'dataParkesNew' filesep  frequency filesep 'results_det_1.mat'], 'output_1090', '-v7.3');
