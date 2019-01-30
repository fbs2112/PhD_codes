import numpy as np
import scipy.signal as sp


class Differentiator:
    def __init__(self, filt_length=3, window=np.hamming(3)):

        if filt_length <= 0:
            raise ValueError('FIR length for differentiator must be positive.')
        if filt_length % 2 != 1:
            raise ValueError('FIR length for differentiator must be odd.')
        
        self.filt_length = filt_length
        self.window = window

        self.idealImpulseResponse = self.get_impulse_response()
        self.truncImpulseResponse = self.get_trunc_impulse_response()
        
    def get_impulse_response(self):

        n = np.arange(1, ((self.filt_length - 1)/2) + 1)
        idealImpulseResponse = (-1**n)/n
        idealImpulseResponse = np.concatenate((-idealImpulseResponse[::-1], idealImpulseResponse))
        idealImpulseResponse = np.insert(idealImpulseResponse, len(idealImpulseResponse)//2, 0)
        return idealImpulseResponse

    def get_trunc_impulse_response(self):
        if len(self.window) != len(self.idealImpulseResponse):
            raise ValueError('Window and impulse response sizes do not match')
        if np.sum(self.window) != 1:
            self.window = self.window/np.sum(self.window)
        return self.idealImpulseResponse*self.window

    def diff_eval(self, x):
        diff = sp.lfilter(self.truncImpulseResponse, 1., x)
        return diff

