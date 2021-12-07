%This function performs the spectrogram evaluation and calls the NMF
%decomposition function.
%Inputs:
%mixtureSignal: input signal. If more than one JNR value is used,
%mixtureSignal must represent a matrix, where each column is composed by
%the received signal for each JNR.
%params: struct with the following mandatory fields:
%        numberOfSources: defines the decomposition rank
%        init: 'random' or 'custom'
%        betaDivergence: either 'kullback-leibler' or 'euclidean'
%        numberOfIterations: max number of iterations
%        tolChange: tolerance parameter that defines convergence
%        tolError: error parameter that defines convergence
% 
%Outputs:
%W2: cell containing the corresponding W matrix for each JNR value
%H2: cell containing the corresponding H matrix for each JNR value
%ReconstructError2: cell containing the corresponding reconstruction error 
%for each JNR value
%dataCell: cell containing the magnitude spectrogram for each JNR value
%f: frequency vector for spectrogram plotting
%t: time vector for spectrogram plotting
%data: cell containing the input signal for each JNR value
%
%This code was created by Felipe Barboza da Silva
% Copyright (c) School of Engineering. Macquarie University - Australia.
% All rights reserved. 2020.
% DISCLAIMER:
%     This code is strictly private, confidential and personal to its recipients
%     and should not be copied, distributed or reproduced in whole or in part,
%     nor passed to any third party.
% ** Do not duplicate or distribute without written permission from the owners. **


function [W2, H2, reconstructError2, dataCell, f, t, data] = nmf_eval_v2(mixtureSignal, params, varargin)

dataCellLength = 1;

if ~isfield(params, 'transform')
    params.transform = true;
end

if ~isfield(params, 'alg')
    params.alg = 'vanilla';
end

if ~isfield(params, 'semi')
    params.semi = false;
end

if ~isfield(params, 'verbose')
    params.verbose = false;
end

if ~isfield(params, 'reassigned')
    params.reassigned = false;
end

if ~isfield(params, 'centered')
    params.centered = true;
end

if ~isfield(params, 'mu')
    params.mu = [0 0];
end

if ~isfield(params, 'transpose')
    params.transpose = false;
end

if nargin > 2 && varargin{1}
    dataCellLength = 2;
end

if ~isfield(params, 'type')
    params.type = 'mag';
end

if ~isfield(params, 'tf')
    params.tf = 'stft';
end

if ~isfield(params, 'betaDivergence')
    params.betaDivergence = 1;
elseif strcmp(params.betaDivergence, 'kullback-leibler')
    params.betaDivergenceAux = 1;
end

dataCell = cell(dataCellLength, length(params.JNRVector));
W2 = cell(dataCellLength, length(params.JNRVector));
H2 = cell(dataCellLength, length(params.JNRVector));
reconstructError2 = cell(dataCellLength, length(params.JNRVector));
data = cell(1, length(params.JNRVector));

for i = 1:size(mixtureSignal, 2)
    
    data{1, i} = mixtureSignal(:,i);
    
    if strcmp(params.tf, 'stft') 
        if params.reassigned
            if params.centered
                [PxxAux, f, t] = spectrogram(data{1, i}, params.window, params.overlap, params.nfft, params.fs, 'centered', params.specType, 'reassigned');
            else
                [PxxAux, f, t] = spectrogram(data{1, i}, params.window, params.overlap, params.nfft, params.fs, params.specType, 'reassigned');
            end
        else
            if params.centered
                [PxxAux, f, t] = spectrogram(data{1, i}, params.window, params.overlap, params.nfft, params.fs, 'centered', params.specType);
            else 
                [PxxAux, f, t] = spectrogram(data{1, i}, params.window, params.overlap, params.nfft, params.fs, params.specType);
            end
        end
    elseif strcmp(params.tf, 'fsst') 
        [PxxAux, f, t] = fsst(data{1, i}, params.fs, params.window);
    else
        error('Time-frequency transform not available');
    end
    
    dataCell{1, i} = PxxAux;
    if params.transpose
        dataCell{1, i} = dataCell{1, i}.';
    end

    if nargin == 3 
        if varargin{1}
            dataCell{2, i} = TK_filtering(PxxAux);
        end
    end
    
    for j = 1:size(dataCell, 1)
        if strcmp(params.type, 'mag')
            inputNMF = abs(dataCell{j, i});
        elseif strcmp(params.type, 'power')
            inputNMF = abs(dataCell{j, i}).^2;
        end
        inputNMF = inputNMF + eps;
        if nargin > 3 && strcmp(varargin{2}, 'filt')
            inputNMF_TPSW = zeros(size(inputNMF));
            for timeIdx = 1:size(inputNMF, 1)
                freqComponent = filter(params.tpsw_coeffs, 1, inputNMF(timeIdx,:));
                inputNMF_TPSW(timeIdx,:) = inputNMF(timeIdx,:);
                inputNMF_TPSW(timeIdx, inputNMF(timeIdx,:) > params.alpha*freqComponent) = freqComponent(inputNMF(timeIdx,:) > params.alpha*freqComponent);
                inputNMF_TPSW(timeIdx,:) = filter(params.tpsw_coeffs, 1, inputNMF_TPSW(timeIdx,:));
            end
            inputNMF = inputNMF - inputNMF_TPSW;
            inputNMF(inputNMF < 0) = eps;
        end
        
        WAux = cell(params.repetitions, 1);
        HAux = cell(params.repetitions, 1);
        reconstructErrorAux = zeros(params.numberOfIterations, params.repetitions);
        
        minReconstructErrorPast = inf;
        
        if params.semi
            wAux = params.W0;
            params.W0 = params.W0(:,(i-1) * params.numberOfSources + 1:i * params.numberOfSources);
        end
        
        for k = 1:params.repetitions
            
            switch params.alg
                case 'vanilla'
                    [W, H, reconstructError] = nmf_vanilla(inputNMF, params);
                case 'vanilla_sparse_H'
                    [W, H, reconstructError] = nmf_vanilla_sparse_H(inputNMF, params);
                case 'vanilla_ort_W'
                    [W, H, reconstructError] = nmf_vanilla_orthogonal_W(inputNMF, params);
                case 'vanilla_ort_H'
                    [W, H, reconstructError] = nmf_vanilla_orthogonal_H(inputNMF, params);
                case 'snmf'
                    [W, H, reconstructError] = nmf_sparse_well_done(inputNMF, params);
                case 'snmf_v2'
                    [W, H, reconstructError] = nmf_sparse_well_done_v2(inputNMF, params);
                case 'vanilla_semi'
                    [W, H, reconstructError] = nmf_vanilla_semi(inputNMF, params);
                case 'snmf_semi'
                    [W, H, reconstructError] = nmf_sparse_well_done_semi(inputNMF, params);
                otherwise
                    error('NMF algorithm not supported');
            end
            WAux{k} = W;
            HAux{k} = H;
            reconstructErrorAux(:,k) = reconstructError;
            idx = find(reconstructErrorAux(:,k), 1, 'last');
            minReconstructErrorPresent = min(reconstructError(1:idx));
            if minReconstructErrorPresent < minReconstructErrorPast
                minReconstructErrorPast = minReconstructErrorPresent;
                minIdx = k;
            end
        end
        W2{j, i} = WAux{minIdx, 1};
        H2{j, i} = HAux{minIdx, 1};
        reconstructError2{j, i} = reconstructErrorAux(:,minIdx);
        
        if params.semi
            params.W0 = wAux;
        end
    end
        
end