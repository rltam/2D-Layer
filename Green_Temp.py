import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
from utility import Parse_vtk_1D
import math
import scipy
from scipy import interpolate
from scipy.interpolate import BSpline, splrep, spalde
import os 
import sys
from Basis_Funcs import Spline_Basis

oldGreens = np.loadtxt('GreenFunc/matrix.txt')
newGreens = np.zeros((len(oldGreens), 31))

offset = np.loadtxt('Diff_Strain_Profile.txt')

for i in range(30):
    newGreens[:,i] = oldGreens[:,i]

newGreens[:,30] = offset

np.savetxt('GreenFunc/matrixOffset.txt', newGreens)