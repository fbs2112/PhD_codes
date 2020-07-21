function [W2, H2, reconstructError2, dataCell, f, t, data] = nmf_eval_v2(mixtureSignal, params, varargin)
%New NMF evaluation function
%Differences from nmf_eval:
    %Input signal may be corrupted by noise

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

if nargin > 2 && varargin{1}
    dataCellLength = 2;
end

dataCell = cell(dataCellLength, length(params.JNRVector));
W2 = cell(dataCellLength, length(params.JNRVector));
H2 = cell(dataCellLength, length(params.JNRVector));
reconstructError2 = cell(dataCellLength, length(params.JNRVector));
data = cell(1, length(params.JNRVector));

for i = 1:size(mixtureSignal, 2)
    
    data{1, i} = mixtureSignal(:,i);
    
    if params.reassigned
        [PxxAux, f, t] = spectrogram(data{1, i}, params.window, params.overlap, params.nfft, params.fs, 'centered', params.specType, 'reassigned');
    else
        [PxxAux, f, t] = spectrogram(data{1, i}, params.window, params.overlap, params.nfft, params.fs, 'centered', params.specType);
    end
    dataCell{1, i} = PxxAux;

    if nargin == 3 
        if varargin{1}
            dataCell{2, i} = TK_filtering(PxxAux);
        end
    end
    
    for j = 1:size(dataCell, 1)
        inputNMF = abs(dataCell{j, i});
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
                case 'snmf'
                    [W, H, reconstructError] = nmf_sparse_well_done(inputNMF, params);
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