import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils import check_array
from sklearn.decomposition.nmf import _initialize_nmf

def normalize_sum(y):
    return y/np.sum(y, 0)


def NMF_denoising_update(data, W_n=None, H_n=None, n_components=2, 
                        init='random',number_iterations=200, 
                        number_init_iterations=20, random_state=1):     


    if init != 'custom':
        W_n, H_n = _initialize_nmf(data, n_components, init=init,
                               random_state=random_state)
        
    interference_data = np.zeros(data.shape)
    interference_data_sum = np.zeros(data.shape)
    
    for i in range(number_iterations):
        W_n *= np.dot((data/(np.dot(W_n, H_n) + interference_data)), H_n.T)
        H_n *= np.dot(W_n.T,(data/(np.dot(W_n, H_n) + interference_data)))

        if i is number_init_iterations:
            aux = (data/np.dot(W_n, H_n))
            interference_data = np.where(aux < 2, 0, aux)
        interference_data *=  data/(np.dot(W_n, H_n) + interference_data) 
        interference_data_sum += interference_data
        data -= interference_data_sum
        W_n = normalize_sum(W_n)

    return W_n, H_n, interference_data_sum


class NMF_denoising(BaseEstimator, TransformerMixin):
    def __init__(self, n_components=2, init='random', number_iterations=200,
                 number_init_iterations=20, random_state=1):
        self.n_components = n_components
        self.init = init
        self.number_iterations = number_iterations
        self.number_init_iterations = number_init_iterations
        self.random_state = random_state

    def fit_transform(self, data, y=None, W_n=None, H_n=None):
        
        data = check_array(data, accept_sparse=('csr', 'csc'), dtype=float)
        W_n, H_n, interference_data_sum = NMF_denoising_update(
            data=data, W_n=W_n, H_n=H_n, n_components=self.n_components, init=self.init,
            number_iterations=self.number_iterations,
            number_init_iterations=self.number_init_iterations,
            random_state=self.random_state)

        self.n_components_ = H_n.shape[0]
        self.components_ = H_n
        self.interference_data = interference_data_sum
        return W_n

    def fit(self, data, y=None, **params):
        self.fit_transform(data, **params)
        return self
