import nmrglue as ng
import pylab
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
 
dic, data = ng.pipe.read("C13_4D.ft1")
sio.savemat('d2o_hcch_tocsy.mat',{'data':data}, format='5')

