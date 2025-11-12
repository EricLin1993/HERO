import nmrglue as ng
import pylab
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
 
dic, data = ng.pipe.read("D2O_Cnoesy.ft1")
sio.savemat('D2O_Cnoesy.mat',{'data':data}, format='5')

