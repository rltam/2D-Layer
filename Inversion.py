import sys
import numpy as np
import netCDF4
import numpy.linalg
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
import pandas as pd
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import sph_harm
import math
import scipy
from scipy.interpolate import BSpline, splrep
from Basis_Funcs import Spline_Basis
import shutil
import os
import datetime
from datetime import datetime
import astropy
from astropy.visualization import hist
from scipy.stats import norm

try:
	n_spline = int(sys.argv[2])
except IndexError:
	n_spline = 25

try:
	deg = int(sys.argv[3])
except IndexError:
	deg = 3

try:
    amp = int(sys.argv[4])
except IndexError:
    amp = 1

try:
    folder = sys.argv[1]
except IndexError:
    folder = 'misc'

try:
    iteration = int(sys.argv[5])
except IndexError:
    iteration = 0

try:
    os.mkdir(f'/home/rltam/Documents/RunOutputFiles/{folder}/Inversion{iteration}')
    dest_folder = f'/home/rltam/Documents/RunOutputFiles/{folder}/Inversion{iteration}'
except FileExistsError:
    dest_folder = f'/home/rltam/Documents/RunOutputFiles/{folder}/Inversion{iteration}'

Edge_Height = 1000
# Edge_Height=np.loadtxt("Edge_Height.txt")
Edge_Length=np.loadtxt("Edge_Length.txt")
true_thick = np.loadtxt(f'/home/rltam/Documents/RunOutputFiles/{folder}/True_Data/Thickness_Profile.txt')
true_thick = true_thick[true_thick[:,0].argsort()]
x = true_thick[:,0]

obs_strain = np.loadtxt(f'/home/rltam/Documents/RunOutputFiles/{folder}/True_Data/Strain_Profile.txt')
try:
    rec_strain = np.loadtxt(f'/home/rltam/Documents/RunOutputFiles/{folder}/{iteration}_Greens({n_spline},{deg})x{amp}/Spline-1/Strain_Profile.txt')
except FileNotFoundError: 
    rec_strain = obs_strain

# G = np.loadtxt(f'GreenFunc/matrix({n_spline},{deg})x{amp}Offset.txt')
G = np.loadtxt(f'GreenFunc/matrix({n_spline},{deg})x{amp}_{iteration}.txt')
G = G[:-1]

if iteration > 0:
    diff_strain = np.zeros(len(obs_strain))
    for i in range(len(diff_strain)):
        diff_strain[i] = obs_strain[i][1] - rec_strain[i][1]
    np.savetxt(f'/home/rltam/Documents/RunOutputFiles/{folder}/Inversion{iteration}/Obs-Rec_Strain.txt', diff_strain)
    
    diff_strain_cut = diff_strain[:-1]

elif iteration == 0:
    diff_strain = np.loadtxt(f'/home/rltam/Documents/RunOutputFiles/{folder}/True_Data/Diff_Strain_Profile.txt')
    diff_strain_cut = diff_strain[:,1][:-1]

basis = Spline_Basis(n_spline, deg)

inverse = scipy.linalg.lstsq(G, diff_strain_cut)
coeffs = inverse[0]

if iteration > 0:
    old_coeffs = np.loadtxt(f'/home/rltam/Documents/RunOutputFiles/{folder}/Inversion{iteration-1}/DerivedCoefficients({n_spline},{deg})x{amp}.txt')
    coeffs = [coeffs[i] + old_coeffs[i] for i in range(len(old_coeffs))]

np.savetxt(f'{dest_folder}/DerivedCoefficients({n_spline},{deg})x{amp}.txt', coeffs)
derived_thick = np.zeros(len(x))

# prev_thick = np.loadtxt(f'/home/rltam/Documents/RunOutputFiles/{folder}/Inversion{iteration-1}/DerivedThickness({n_spline},{deg})x{amp}.txt')

for i in range(len(x)):
    val = 0
    for j in range(len(basis)):
        base = basis[j]
        if j == len(basis) - 1:
            if x[i] >= 49999:
                val += Edge_Height * amp/100 *coeffs[j]
            else:
                val += Edge_Height * amp/100*coeffs[j]*max(base(x[i]), 0)
        else:
            val += Edge_Height * amp/100 *coeffs[j]*max(base(x[i]), 0)
    # val += Edge_Height * amp/100 * coeffs[-1]

    # Val = prev_thick[:,1][i] + 2*val
    Val = 1000 + 2*val
    derived_thick[i] = Val

derived_prof = np.zeros((len(x), 2))
derived_prof[:,0] = x
derived_prof[:,1] = derived_thick

np.savetxt(f'{dest_folder}/DerivedThickness({n_spline},{deg})x{amp}.txt', derived_prof)

fig, ax = plt.subplots(figsize = (10,7))
ax.plot(x, true_thick[:,1], c='k', label = "True")
for i in range(iteration):
    old_thick = np.loadtxt(f'/home/rltam/Documents/RunOutputFiles/{folder}/Inversion{i}/DerivedThickness({n_spline},{deg})x{amp}.txt')
    ax.plot(x, [old_thick[:,1][i] for i in range(len(x))], label=f'Derived (Iteration {i})')
ax.plot(x, [derived_thick[i] for i in range(len(x))], label=f"Derived (Iteration {iteration})")
# ax.scatter(x, derived_thick, c='r', s=10, label="Derived")
# ax.scatter(x, true_thick[:,1], c='b', s=10, label = "True")
ax.legend()
ax.set_xlabel("Distance Along Block (m)")
ax.set_ylabel("Thickness Profile from (m)")
ax.set_title(f'True vs Derived Thickness Profiles')
plt.savefig(f'{dest_folder}/ProfilePlot({n_spline},{deg})x{amp}.png')
plt.close()

resids = np.zeros((len(x),5))
resids[:,0] = x
resids[:,1] = [derived_thick[i] - true_thick[:,1][i] for i in range(len(x))]
resids[:,2] = [resids[:,1][i]/true_thick[:,1][i] for i in range(len(x))]
resids[:,3] = [(resids[:,1][i])**2 for i in range(len(x))]
resids[:,4] = [(resids[:,2][i])**2 for i in range(len(x))]

np.savetxt(f'{dest_folder}/Residuals({n_spline},{deg})x{amp}.txt', resids)

ave = np.zeros(4)
for i in range(4):
    ave[i] = np.mean(resids[:,i+1])

np.savetxt(f'{dest_folder}/AveResiduals({n_spline},{deg})x{amp}.txt', ave)

shutil.copy(f'{dest_folder}/DerivedCoefficients({n_spline},{deg})x{amp}.txt', 'DerivedCoefficients.txt')
