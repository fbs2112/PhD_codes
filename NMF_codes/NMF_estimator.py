from sklearn.decomposition import NMF
from sklearn.decomposition.nmf import _beta_divergence
 
 
class NMF_estimator(NMF):
     
    # def __init__(self, n_components=None, init=None, solver='cd',
    #              beta_loss='frobenius', tol=1e-4, max_iter=200,
    #              random_state=None, alpha=0., l1_ratio=0., verbose=0,
    #              shuffle=False):
    #     super().__init__(n_components, init, solver,
    #              beta_loss, tol, max_iter,
    #              random_state, alpha, l1_ratio, verbose,
    #              shuffle)  
 
    def fit_transform(self, X, y=None, W=None, H=None):
        W = super().fit_transform(X, y, W, H)
        self.W_ = W
 
        return W
     
 
    def score(self, X, y=None):
        H = self.components_
        W = self.transform(X)
        return -_beta_divergence(X, W, H, self.beta_loss, square_root=True)
   

