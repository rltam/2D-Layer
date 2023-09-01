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


strain_profile = np.loadtxt('output/Strain_Profile.txt')
base = np.loadtxt('GreenFunc/x_strain_unperturbed.txt')
strain = strain_profile

frac_strain = [strain[i] - base[i] for i in range(len(base))]

G = np.loadtxt('GreenFunc/matrix.txt')

## cut out edges:
G_cut = G[:-1]
frac_strain_cut = frac_strain[:-1]

inverse = scipy.sparse.linalg.lsqr(G_cut, frac_strain_cut, damp = 5e-8, show=False)

coeffs = inverse[0]
print(coeffs)
basis = Spline_Basis(30,3)

Edge_Height=np.loadtxt("Edge_Height.txt")
Edge_Length=np.loadtxt("Edge_Length.txt")

true_thick = np.loadtxt('output/Thickness_Profile.txt')

true_thick = true_thick[true_thick[:,0].argsort()]
x = true_thick[:,0]
# x = np.linspace(-Edge_Length/2, Edge_Length/2, 2000)

derived_prof = []

for i in range(len(x)):
    val = 0
    for j in range(len(basis)):
        base = basis[j]
        if j == len(basis) - 1:
            if x[i] >= 49999:
                val += 0.9998800047999359 * Edge_Height * 0.01*coeffs[j]
            else:
                val += Edge_Height * 0.01*coeffs[j]*max(base(x[i]), 0)
        else:
            val += Edge_Height * 0.01*coeffs[j]*max(base(x[i]), 0)
    
    derived_prof.append(val)


fig, ax1 = plt.subplots()
# ax2 = ax1.twinx()

# ax1.scatter(true_thick[:,0], [true_thick[:,1][i]-1000 for i in range(len(true_thick[:,1]))], s=5, c = 'r', label = 'True Thickness Profile')
ax1.plot(true_thick[:,0], [true_thick[:,1][i]-1000 for i in range(len(true_thick[:,1]))], c = 'r', label = 'True Thickness Profile')
ax1.plot(x, derived_prof, label = 'Derived Profile')
# for i in range(len(basis)):
        # plt.plot(x, [100*coeffs[i]*basis[i](x[j]) for j in range(len(x))])
ax1.legend(loc="upper left")
# ax2.legend(loc="upper right")
plt.show()
    

## residuals - chi squared
diffs = np.zeros(len(x))
diffs_squared = np.zeros(len(x))
## should we also normalize by the scale of deviations in thickness from uniform model???

for i in range(len(x)):
    diffs[i] = derived_prof[i] - (true_thick[:,1][i] - 1000)
    diffs_squared[i] = (derived_prof[i] - (true_thick[:,1][i] - 1000))**2

sum_RS = np.sum(diffs_squared) / len(x)

