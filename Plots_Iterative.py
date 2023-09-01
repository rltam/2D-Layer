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

Edge_Height=np.loadtxt("Edge_Height.txt")
Edge_Length=np.loadtxt("Edge_Length.txt")

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
    iteration = int(sys.argv[5]) + 1
except IndexError:
    iteration = 11

dest_folder = f'/home/rltam/Documents/RunOutputFiles/{folder}'
left = -Edge_Length/2
right = Edge_Length/2

basis = Spline_Basis(n_spline, deg)

real_coeffs = np.loadtxt(f'{dest_folder}/True_Data/Coefficients.txt')
true_strain = np.loadtxt(f'{dest_folder}/True_Data/Strain_Profile.txt')

x = np.linspace(left, right, 2000)
true_thick = np.zeros(len(x))

mean_thick_list = []

for i in range(len(x)):
    val = 0
    for j in range(len(basis)):
        base = basis[j]
        if j == len(basis) - 1:
            if x[i] >= 49999:
                val += Edge_Height * amp/100 *real_coeffs[j]
            else:
                val += Edge_Height * amp/100*real_coeffs[j]*max(base(x[i]), 0)
        else:
            val += Edge_Height * amp/100 *real_coeffs[j]*max(base(x[i]), 0)

    Val = Edge_Height + 2*val
    true_thick[i] = Val

mean_thick_list.append(np.mean(true_thick))
diff_list = []
thickness_resid_list = []
strain_resid_list = []

# fig, ax = plt.subplots(figsize=(15,10))
# ax.plot(x, true_thick, c='k', label='True')

for n in range(iteration):
    coeffs = np.loadtxt(f'{dest_folder}/Inversion{n}/DerivedCoefficients({n_spline},{deg})x{amp}.txt')
    derived_thick = np.zeros(len(x))
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

        Val = 1000 + 2*val
        derived_thick[i] = Val
    # mean_thick_list.append(np.mean(derived_thick))
    # diff_list.append(mean_thick_list[0] - np.mean(derived_thick))
    derived_profile = np.zeros((len(x), 2))
    derived_profile[:,0] = x
    derived_profile[:,1] = derived_thick
    np.savetxt(f'{dest_folder}/Inversion{n}/Derived_Thickness_Profile.txt', derived_profile)

    resids = np.zeros((len(x),5))
    resids[:,0] = x
    resids[:,1] = [np.abs(derived_thick[i] - true_thick[i]) for i in range(len(x))]
    resids[:,2] = [np.abs(derived_thick[i] - true_thick[i])/true_thick[i] for i in range(len(x))]
    resids[:,3] = [(resids[:,1][i])**2 for i in range(len(x))]
    resids[:,4] = [(resids[:,2][i])**2 for i in range(len(x))]
    np.savetxt(f'{dest_folder}/Inversion{n}/Thickness_Residuals_{n}.txt', resids)
    thickness_resid_list.append(np.sqrt(np.mean(resids[:,4])))

    if n < iteration-1:
        strain_diff = np.loadtxt(f'{dest_folder}/Inversion{n+1}/Obs-Rec_Strain.txt')

    else:
        strain_diff = np.loadtxt(f'{dest_folder}/Last_Iteration/Obs-Rec_Strain.txt')
    
    
    resids_st = np.zeros((len(strain_diff), 4))
    resids_st[:,0] = [np.abs(strain_diff[i]) for i in range(len(strain_diff))]
    resids_st[:,1] = [np.abs(strain_diff[i])/true_strain[i] for i in range(len(strain_diff))]
    resids_st[:,2] = [(resids_st[:,0][i])**2 for i in range(len(strain_diff))]
    resids_st[:,3] = [(resids_st[:,1][i])**2 for i in range(len(strain_diff))]
    resids_st = resids_st[:-1]
    np.savetxt(f'{dest_folder}/Inversion{n}/Strain_Residuals_{n}.txt', resids)
    strain_resid_list.append(np.sqrt(np.mean(resids_st[:,3])))


    fig, ax = plt.subplots(figsize=(10,7))
    ax.plot(x, true_thick, c='k', label='True')
    ax.plot(x, derived_thick[i], label = 'Derived')
    ax.set_xlabel('Distance Along Block (m)')
    ax.set_ylabel('Thickness (m)')
    ax.set_title(f'True vs Derived Thickness Profiles from Iteration {n}')
    ax.legend()
    plt.savefig(f'{dest_folder}/Thickness_Profile_{n}.png')
    # shutil.copy(f'Thickness_Profile_{n}.png', f'Inversion{n}/Shifted_Thickness_Profile.png')
    plt.close()

    fig, ax = plt.subplots(2, 1, sharex = True, figsize=(15,10))
    ax[0].plot(x, true_thick, c='k', label='True')
    ax[0].plot(x, derived_thick[i], label = 'Derived')
    ax[0].set_xlabel('Distance Along Block (m)')
    ax[0].set_ylabel('Thickness (m)')
    ax[0].set_title(f'True vs Derived Thickness Profiles from Iteration {n}')
    ax[0].legend()

    ax2 = ax[1].twinx()
    ax[1].plot(x, resids[:,2], c = 'r', label=f'Thickness')    
    y = np.linspace(-5e4, 5e4, len(strain_diff))
    ax2.plot(y, resids_st[:,1], c='b', label=f'Strain')
    ax[1].set_xlabel('Distance Along Block (m)')
    ax[1].set_ylabel('Thickness Fractional Error')
    ax2.set_ylabel('Strain Fractional Error')
    ax[1].legend(loc="upper left")
    ax2.legend(loc="upper right")
    ax[1].set_title('Fractional Error of Recovered Profile')
    plt.savefig(f'{dest_folder}/Inversion{n}/Profiles+Residuals.png')
    plt.close()

# ax.plot(x, true_thick, c='k')
# ax.set_xlabel('Distance Along Block (m)')
# ax.set_ylabel('Thickness (m)')
# ax.set_title(f'True vs Derived Thickness Profiles')
# ax.legend()
# #plt.show()
# plt.savefig('Thickness_Profiles.png')
# plt.close()

# last_strain = np.loadtxt(f'{folder}/Last_Iteration/Obs-Rec_Strain.txt')
# last_strain = [last_strain[i]**2 for i in range(len(last_strain))]
# strain_resid_list.append(np.mean(last_strain))

#recent_strain = np.loadtxt('24_Greens(25,3)x1/Spline-1/Strain_Profile.txt')
#diff_last_strain = [(true_strain[:,1][i] - recent_strain[:,1][i])**2 for i in range(len(true_strain))]
#strain_resid_list.append(np.mean(diff_last_strain))

iteration = np.linspace(0,iteration-1,iteration)
fig, ax1 = plt.subplots(figsize=(10,7))
ax2 = ax1.twinx()
ax1.plot(iteration, thickness_resid_list, c='r', label='Thickness')
ax2.plot(iteration, strain_resid_list, c='b', label='Strain')
ax1.set_xlabel('Inversion Iteration Number')
ax1.set_ylabel('Thickness Fractional Error RMS')
ax2.set_ylabel('Strain Fractional Error RMS')
ax1.legend(loc="upper left")
ax2.legend(loc="upper right")
ax1.set_title('RMS Residuals by Iteration')
plt.savefig(f'{dest_folder}/Residuals.png')
plt.close()