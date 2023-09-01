import sys
import numpy as np
import netCDF4
import numpy as np
import numpy.linalg
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
import scipy
from scipy import interpolate
from scipy.interpolate import BSpline, splrep
from Basis_Funcs import Spline_Basis

n_spline = int(sys.argv[1])

def get_coefficients(n):
    coeffs = np.random.uniform(-1, 1, n)
    return coeffs

coeffs = get_coefficients(n_spline)
# print(coeffs)
np.savetxt('Coefficients.txt', coeffs)
