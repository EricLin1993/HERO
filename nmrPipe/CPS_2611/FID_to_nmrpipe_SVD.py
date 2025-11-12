import nmrglue as ng
import pylab
import numpy as np
import scipy.io as sio


dic, data = ng.pipe.read("csr4_n15_noesy.ft1")
M = sio.loadmat('Nnoesy_recon_SVD.mat')
Data = M['FID_Rec']
Data = np.array(Data,dtype='float32')
ng.pipe.write("Nnoesy_recon_SVD.ft3", dic, Data, overwrite=True)
