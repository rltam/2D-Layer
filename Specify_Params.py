



# *********************************************************************
# FUNCTION TO DEFINE PARAMETERS AT START OF SIMULATION
# *********************************************************************


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
import pandas as pd
import pyshtools as pysh
from pyshtools import constants
from pyshtools import gravmag
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import sph_harm
import math
import pickle 



Edge_Length = 1e5 #m


Edge_Height = 1e3 #m



Edge_Length=np.array([Edge_Length])
np.savetxt("Edge_Length.txt",Edge_Length,fmt ="%d")
tstr=str(Edge_Length)
print("Edge_Length:")
print(tstr)


Edge_Height=np.array([Edge_Height])
np.savetxt("Edge_Height.txt",Edge_Height,fmt ="%d")
tstr=str(Edge_Height)
print("Edge_Height:")
print(tstr)



