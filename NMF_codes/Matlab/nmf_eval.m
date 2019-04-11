
function [W2, H2, reconstructError2, dataCell, f, t, data] = nmf_eval(mixtureSignal, params, varargin)

dataCellLength = 1;

if nargin > 2 && varargin{1}
    dataCellLength = 2;
end

dataCell = cell(dataCellLength, length(params.JNRVector));
W2 = cell(dataCellLength, length(params.JNRVector));
H2 = cell(dataCellLength, length(params.JNRVector));
reconstructError2 = cell(dataCellLength, length(params.JNRVector));
data = cell(1, length(params.JNRVector));

powMixture = pow_eval(mixtureSignal);
signal1Length = length(mixtureSignal);

noise = randn(signal1Length, 1);

if ~isreal(mixtureSignal)
    noise = randn(signal1Length, 1) + 1j*randn(signal1Length, 1);
end

powNoise = pow_eval(noise);

for i = 1:length(params.JNRVector)
    
    powAux = powMixture/db2pow(params.JNRVector(i));
    noise2 = noise*sqrt(powAux/powNoise);
    data{1, i} = mixtureSignal + noise2;
    
    if isreal(data{1, i})
        [PxxAux, f, t] = spectrogram(data{1, i}, hann(params.nperseg), params.overlap, params.nfft, params.fs, 'power');
    else
        [PxxAux, f, t] = spectrogram(data{1, i}, hann(params.nperseg), params.overlap, params.nfft, params.fs, 'centered', 'power');
    end    
    
    dataCell{1, i} = PxxAux;

    if nargin > 2 
        if varargin{1}
            dataCell{2, i} = TK_filtering(PxxAux);
        end
    end
    
    for j = 1:size(dataCell, 1)
        inputNMF = abs(dataCell{j, i}).^2;
        WAux = cell(params.repetitions, 1);
        HAux = cell(params.repetitions, 1);
        reconstructErrorAux = zeros(params.numberOfIterations, params.repetitions);
        
        minReconstructErrorPast = inf;
        for k = 1:params.repetitions
            [W, H, reconstructError] = nmf_py3(inputNMF, params);
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
        
    end
        
end