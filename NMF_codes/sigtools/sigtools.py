import numpy as np
from scipy.signal import butter, lfilter

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y
    

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y


def db2pow(y):
    return 10.**(y/10.)

def pow_eval(y):
    return (np.linalg.norm(y)**2)/len(y)

def TK_filtering(x):
    
    def TK_eval(x):
        b = x**2
        aux = np.concatenate((np.zeros(1), x[:-1]))
        aux2 = np.concatenate((x[1:], np.zeros(1)))
        a = b - aux2*aux
        return a
    
    isreal = np.isreal(x).all()
    
    if isreal:
        return TK_eval(x)
    else:
        return TK_eval(np.real(x)) + TK_eval(np.imag(x))