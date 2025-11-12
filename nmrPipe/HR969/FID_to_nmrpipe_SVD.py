import nmrglue as ng
import pylab
import numpy as np
import scipy.io as sio


dic, data = ng.pipe.read("Hnoesy_4D.ft1")
M = sio.loadmat('Hnoesy_4D_recon_SVD.mat')
Data = M['FID_Rec']
Data = np.array(Data,dtype='float32')
ng.pipe.write("Hnoesy_4D_recon_SVD.ft3", dic, Data, overwrite=True)
