import nmrglue as ng
import pylab
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
 
dic, data = ng.pipe.read("projYZ.dat")
sio.savemat('projYZ_SVD.mat',{'data':data}, format='5')

