

##Plotting Results....

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
from utility import Parse_vtk_1D

Thick_1=np.loadtxt('Thickness_Profile_n0s.txt')

Thick_2=np.loadtxt('Thickness_Profile_true.txt')

plt.scatter(Thick_1[:,0],Thick_1[:,1])
plt.scatter(Thick_2[:,0],Thick_2[:,1])
plt.show()







