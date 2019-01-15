import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from scipy.signal import stft, istft

class stft_fun(BaseEstimator, TransformerMixin):
    
    def __init__(self, fs=None, nperseg=None, nfft=None):  
        self.fs = fs
        self.nperseg = nperseg
        self.nfft = nfft


    def fit(self, X=None, y=None):
        return self

    def transform(self, data):
        self.isreal = np.isreal(data).any()
        f, t, Zxx = stft(data, fs=self.fs, nperseg=self.nperseg, nfft=self.nfft, boundary=None)
        self.freq = f
        self.time = t
        self.Zxx = Zxx
        return np.abs(Zxx)

    def inverse_transform(self, W, H):
        numberOfSources = W.shape[1]
        self.phase = np.angle(self.Zxx)
        dataHat = []
        for i in range(numberOfSources):
            x = np.outer(W[:,i],H[i,:])
            x = x*np.exp(1j*self.phase)
            _, y = istft(x, fs=self.fs, nperseg=self.nperseg, nfft=self.nfft, input_onesided=self.isreal, boundary=False)
            dataHat.append(y)
        return dataHat
   

   