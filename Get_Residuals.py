import sys
import numpy as np
import netCDF4
import numpy as np
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
	iter_num = sys.argv[1]
except IndexError:
	iter_num = ''

# saving files etc
txtfiles = ['Displacement_Profile', 'Strain_Profile', 'Thickness_Profile', 'Diff_Strain_Profile']
vtkfiles = ['UpperSurf', 'LeftSurf', 'RightSurf', 'All_Nodes']

dt_string = sys.argv[1]

os.mkdir(f'/home/rltam/Documents/RunOutputFiles/Residual_Analysis/run{iter_num}')


dest_folder = f'/home/rltam/Documents/RunOutputFiles/Residual_Analysis/run{iter_num}'

for file in txtfiles:
    shutil.copy(f'output/{file}.txt', f'{dest_folder}/{file}.txt')

for file in vtkfiles:
    shutil.copy(f'output/{file}_t0000000.vtk', f'{dest_folder}/{file}.vtk')

shutil.copy('Coefficients.txt', f'{dest_folder}/TrueCoefficients.txt')

# everything else

Edge_Height=np.loadtxt("Edge_Height.txt")
Edge_Length=np.loadtxt("Edge_Length.txt")
true_thick = np.loadtxt('output/Thickness_Profile.txt')
# input_coeffs = np.loadtxt('Coefficients.txt')
frac_strain = np.loadtxt('Diff_Strain_Profile.txt')

G = np.loadtxt('GreenFunc/matrix50.txt')

## cut out edges:
G_cut = G[:-1]
frac_strain_cut = frac_strain[:-1]

basis = Spline_Basis(50,3)

true_thick = true_thick[true_thick[:,0].argsort()]
x = true_thick[:,0]
true_thickness = [true_thick[:,1][i] - 1000 for i in range(len(true_thick[:,1]))]

damping_range = np.logspace(-9, -6, 4)
damping_range = [1e-8, 2.154e-8, 4.641e-8, 1e-7, 2.154e-7, 4.641e-7, 1e-6, 2.154e-6, 4.641e-6, 1e-5]
resid_list = []

print('\n')

for damping in damping_range:
    print(f'Run {iter_num}, damp = {damping}')
    inverse = scipy.sparse.linalg.lsqr(G_cut, frac_strain_cut, damp = damping, show=False)
    coeffs = inverse[0]
    
    np.savetxt(f'{dest_folder}/DerivedCoefficients_{damping}', coeffs)
    derived_prof = []

    for i in range(len(x)):
        val = 0
        for j in range(len(basis)):
            base = basis[j]
            if j == len(basis) - 1:
                if x[i] >= 49999:
                    val += 1 * Edge_Height * 0.01*coeffs[j]
                else:
                    val += Edge_Height * 0.01*coeffs[j]*max(base(x[i]), 0)
            else:
                val += Edge_Height * 0.01*coeffs[j]*max(base(x[i]), 0)
        
        derived_prof.append(val)

    ## residuals 

    diffs = np.zeros(len(x))

    for i in range(len(x)):
        diffs[i] = (derived_prof[i] - true_thickness[i]) / (true_thick[:,1][i])

    diffs_squared = [diffs[i]**2 for i in range(len(diffs))]

    sum_RS = np.mean(diffs_squared)

    np.savetxt(f'{dest_folder}/Residuals_{damping}.txt',diffs)
    np.savetxt(f'{dest_folder}/Residuals2_{damping}.txt',diffs_squared)
    resid_list.append(sum_RS)

    fix, ax = plt.subplots(figsize=(10, 7))

    ax.plot(x, true_thickness, c = 'k', label = 'True Profile')
    ax.plot(x, derived_prof, c = 'r', label = 'Derived Profile')
    ax.legend()
    ax.set_xlabel("Distance Along Block (m)")
    ax.set_ylabel("Thickness Variation from (m)")
    ax.set_title(f'True vs Derived Profile for Run {iter_num}, damp = {damping}')

    plt.savefig(f'{dest_folder}/ThicknessPlot_{damping}.png')

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.scatter(x, diffs, s = 5, c = 'darkcyan', label='Residuals')
    ax.grid()
    ax.set_xlabel("Distance Along Block (m)")
    ax.set_ylabel('Residual')
    ax.set_title(f'Residuals for Run {iter_num}, damp = {damping}, Mean Squared Residual = {sum_RS}')
    plt.savefig(f'{dest_folder}/ResidudalPlot_{damping}.png')


    # make plot
    # fig, ax = plt.subplots(2,2, figsize=(20, 15))

    # ax[0][0].plot(x, true_thickness[:,1], c = 'k', label = 'True Thickness Profile')
    # ax[0][0].plot(x, derived_prof, c = 'r', label = 'Derived Profile')
    # ax[0][0].legend()
    # ax[0][0].set_xlabel("Distance Along Block (m)")
    # ax[0][0].set_ylabel("Thickness Variation (m)")

    # ax[0][1].scatter(x, diffs, s = 8, c = 'darkcyan', label='Residuals')
    # # ax[0][1].scatter(x, diffs_squared, s = 8, c = 'green', label='Squared Residuals')
    # ax[0][1].grid()
    # ax[0][1].set_xlabel("Distance Along Block (m)")
    # ax[0][1].set_ylabel('Residual') 
    # ax[0][1].legend()

    # line = np.linspace(1,30,30)

    # # hist(diffs, bins = 40, ax = ax[1][0], color = 'mediumturquoise', label='Residuals')
    # hist(diffs_squared, bins = 40, ax = ax[1][1], color = 'darkcyan', label='Squared Residuals')
    # # xmin, xmax = ax[1].set_xlim()
    # # x_pdf = np.linspace(xmin, xmax, 100)
    # # p = norm.pdf(x_pdf, mu, std)

    # # ax[1].plot(x_pdf, p, label=r'Gaussian Fit to Squared Residuals Histogram')
    # ax[1][0].scatter(line, coeff_resids, c = 'green', label='Coefficient Residuals')
    # ax[1][0].grid()
    # ax[1][0].legend()
    # ax[1][0].set_xlabel('Coefficient Number')
    # ax[1][0].set_ylabel('Residual') 
    # ax[1][0].legend()
    # ax[1][1].set_xlabel('Value')
    # ax[1][1].set_ylabel('Number of Points') 
    # ax[1][1].legend()

    # fig.suptitle(f'Profile and Residual Histogram for Inverse Method, damp = {damping}, Mean Squared Residual = {str(np.mean(diffs_squared))[:9]}')

    # # plt.show()
    # plt.savefig(f'{dest_folder}/plots_{damping}.png')



sum_data = np.zeros((len(damping_range), 2))

sum_data[:,0] = damping_range
sum_data[:,1] = resid_list

np.savetxt(f'{dest_folder}/AveResiduals.txt', sum_data)

fig, ax1 = plt.subplots()

ax1.scatter(damping_range, resid_list, label='Mean Squared Residual')
ax1.set_xlabel('Damping Factor')
ax1.set_ylabel('Mean Squared Residual')
ax1.set_xscale('log')
ax1.set_yscale('log')


# ax.set_title('Mean Squared Residual vs Damping Factor')
plt.savefig(f'{dest_folder}/MRSvsDamp.png')

# fig, ax = plt.subplots()

# ax.scatter(damping_range, r_squared_list)
# ax.set_xlabel('Damping Factor')
# ax.set_ylabel(r'$r^2$')
# ax.set_title(r'$r^2$ vs Damping Factor')
# plt.savefig(f'{dest_folder}/r2vsDamp.png')
