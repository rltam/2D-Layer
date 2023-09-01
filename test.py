import sys
import numpy as np
import netCDF4
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
import scipy
from scipy import interpolate
from scipy.interpolate import BSpline, splrep
from Basis_Funcs import Spline_Basis

folder = sys.argv[1]
n_splines = int(sys.argv[2])
deg = int(sys.argv[3])

basis = Spline_Basis(n_splines, deg)
print(f'{folder}')

# fig, ax = plt.subplots(2, 2, sharex=True)

# x = np.linspace(-50000, 50000, 100)

# for i in range(len(basis)):
#     base = basis[i]
#     plt.plot(x, [max(0,base(x[i])) for i in range(len(x))])

# plt.show()