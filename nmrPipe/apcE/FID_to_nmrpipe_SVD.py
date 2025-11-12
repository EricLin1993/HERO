import nmrglue as ng
import pylab
import numpy as np
import scipy.io as sio


dic, data = ng.pipe.read("d2o_hcch_tocsy.ft1")
M = sio.loadmat('tocsy_recon_SVD.mat')
Data = M['FID_Rec']
Data = np.array(Data,dtype='float32')
ng.pipe.write("tocsy_recon_SVD.ft3", dic, Data, overwrite=True)
