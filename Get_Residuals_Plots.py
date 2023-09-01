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


# saving files etc
# txtfiles = ['Displacement_Profile', 'Strain_Profile', 'Thickness_Profile', 'Diff_Strain_Profile']
# vtkfiles = ['UpperSurf', 'LeftSurf', 'RightSurf', 'All_Nodes']

# dt_string = folder = sys.argv[1]

# os.mkdir(f'/home/rltam/Documents/RunOutputFiles/{dt_string}')
# os.mkdir(f'/home/rltam/Documents/RunOutputFiles/{dt_string}/RunData')
# os.mkdir(f'/home/rltam/Documents/RunOutputFiles/{dt_string}/InverseOutputs')



# dest_folder = f'/home/rltam/Documents/RunOutputFiles/{dt_string}'

# for file in txtfiles:
#     txtSplit = file.split('_')
#     word = txtSplit[0]
#     shutil.copy(f'output/{file}.txt', f'{dest_folder}/RunData/{word}_{dt_string}.txt')

# for file in vtkfiles:
#     txtSplit = file.split('_')
#     word = txtSplit[0]
#     shutil.copy(f'output/{file}_t0000000.vtk', f'{dest_folder}/RunData/{word}_{dt_string}.vtk')

# shutil.copy('Coefficients.txt', f'{dest_folder}/RunData/Coefficients_{dt_string}.txt')

# everything else

strain = np.loadtxt('output/Strain_Profile.txt')
base = np.loadtxt('GreenFunc/x_strain_unperturbed.txt')
Edge_Height=np.loadtxt("Edge_Height.txt")
Edge_Length=np.loadtxt("Edge_Length.txt")
true_thick = np.loadtxt('output/Thickness_Profile.txt')
input_coeffs = np.loadtxt('Coefficients.txt')


frac_strain = [strain[i] - base[i] for i in range(len(base))]

G = np.loadtxt('GreenFunc/matrixOffset.txt')

## cut out edges:
G_cut = G[:-1]
frac_strain_cut = frac_strain[:-1]

G_30 = np.delete(G_cut, -1, axis = 1)

G_offset = np.zeros((1000, 2))
G_offset[:,0] = G_cut[:,-1]


basis = Spline_Basis(50,3)

true_thick = true_thick[true_thick[:,0].argsort()]
x = true_thick[:,0]

damping_range = np.linspace(5e-8, 5e-7, 10)
resid_list = []
r_squared_list = []
coeff_resid_list = []

for damping in damping_range:
    print(f'Now using damping factor {damping}')
    inverse = scipy.sparse.linalg.lsqr(G_cut, frac_strain_cut, damp = damping, show=False)
    coeffs = inverse[0]

    # inverse1 = scipy.sparse.linalg.lsqr(G_30, frac_strain_cut, damp = damping, show = False)
    # coeffs1 = inverse1[0]

    # strain_no_offset = np.dot(G_30, coeffs1)

    # strain_offset = frac_strain_cut - strain_no_offset

    # inverse2 = np.linalg.lstsq(G_offset, strain_offset, rcond=None)
    # coeffs2 = inverse2[0]
    
    # coeffs = np.concatenate((coeffs1, coeffs2))
    
    np.savetxt(f'{dest_folder}/InverseOutputs/Inverse_Coeffs_{damping}', coeffs)
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
        
        val += Edge_Height * 0.01 * coeffs[-1]

        derived_prof.append(val)

    ## residuals 

    diffs = np.zeros(len(x))
    diffs_squared = np.zeros(len(x))
    ## should we also normalize by the scale of deviations in thickness from uniform model???

    true_thickness = [true_thick[:,1][i] - 1000 for i in range(len(true_thick[:,1]))]

    for i in range(len(x)):
        diffs[i] = derived_prof[i] - true_thickness[i]
        diffs_squared[i] = (derived_prof[i] -  true_thickness[i])**2
        # diffs[i] = (derived_prof[i] - (true_thick[:,1][i] - 1000)) / (true_thick[:,1][i] - 1000)
        # diffs_squared[i] = ((derived_prof[i] - (true_thick[:,1][i] - 1000)) / (true_thick[:,1][i] - 1000))**2

    sum_RS = np.sum(diffs_squared) / len(x)
    np.savetxt(f'{dest_folder}/InverseOutputs/Residuals_{damping}.txt',diffs)
    np.savetxt(f'{dest_folder}/InverseOutputs/Residuals2_{damping}.txt',diffs_squared)
    resid_list.append(sum_RS)

    # r_squared = 1 - np.sum(diffs_squared) / np.sum([(true_thickness[i] - np.mean(true_thickness))**2 for i in range(len(true_thickness))])
    # r_squared_list.append(r_squared)

    coeff_resids = [coeffs[i] - input_coeffs[i] + coeffs[len(input_coeffs)] for i in range(len(input_coeffs))]
    coeff_resids_squared = [coeff_resids[i]**2 for i in range(len(input_coeffs))]
    np.savetxt(f'{dest_folder}/InverseOutputs/CoeffResiduals_{damping}.txt',diffs)
    np.savetxt(f'{dest_folder}/InverseOutputs/CoeffResiduals2_{damping}.txt',diffs_squared)

    sum_CRS = np.mean(coeff_resids_squared)
    coeff_resid_list.append(sum_CRS)

    # make plot
    fig, ax = plt.subplots(2,2, figsize=(20, 15))

    # ax1.scatter(true_thick[:,0], [true_thick[:,1][i]-1000 for i in range(len(true_thick[:,1]))], s=5, c = 'r', label = 'True Thickness Profile')
    ax[0][0].plot(true_thick[:,0], [true_thick[:,1][i]-1000 for i in range(len(true_thick[:,1]))], c = 'r', label = 'True Thickness Profile')
    ax[0][0].plot(x, derived_prof, c = 'k', label = 'Derived Profile')
    ax[0][0].legend()
    ax[0][0].set_xlabel("Distance Along Block (m)")
    ax[0][0].set_ylabel("Thickness Variation (m)")

    ax[0][1].scatter(x, diffs, s = 8, c = 'darkcyan', label='Residuals')
    # ax[0][1].scatter(x, diffs_squared, s = 8, c = 'green', label='Squared Residuals')
    ax[0][1].grid()
    ax[0][1].set_xlabel("Distance Along Block (m)")
    ax[0][1].set_ylabel('Residual') 
    ax[0][1].legend()

    line = np.linspace(1,30,30)

    # hist(diffs, bins = 40, ax = ax[1][0], color = 'mediumturquoise', label='Residuals')
    hist(diffs_squared, bins = 40, ax = ax[1][1], color = 'darkcyan', label='Squared Residuals')
    # xmin, xmax = ax[1].set_xlim()
    # x_pdf = np.linspace(xmin, xmax, 100)
    # p = norm.pdf(x_pdf, mu, std)

    # ax[1].plot(x_pdf, p, label=r'Gaussian Fit to Squared Residuals Histogram')
    ax[1][0].scatter(line, coeff_resids, c = 'green', label='Coefficient Residuals')
    ax[1][0].grid()
    ax[1][0].legend()
    ax[1][0].set_xlabel('Coefficient Number')
    ax[1][0].set_ylabel('Residual') 
    ax[1][0].legend()
    ax[1][1].set_xlabel('Value')
    ax[1][1].set_ylabel('Number of Points') 
    ax[1][1].legend()

    fig.suptitle(f'Profile and Residual Histogram for Inverse Method, damp = {damping}, Mean Squared Residual = {str(np.mean(diffs_squared))[:9]}')

    # plt.show()
    plt.savefig(f'{dest_folder}/plots_{damping}.png')



sum_data = np.zeros((len(damping_range), 2))
coeff_sum_data = np.zeros((len(damping_range), 2))

sum_data[:,0] = damping_range
sum_data[:,1] = resid_list
coeff_sum_data[:,0] = damping_range
coeff_sum_data[:,1] = coeff_resid_list

np.savetxt(f'{dest_folder}/InverseOutputs/AveResiduals.txt', sum_data)
np.savetxt(f'{dest_folder}/InverseOutputs/AveCoeffResiduals.txt', coeff_sum_data)

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()  

ax1.scatter(damping_range, resid_list, label='Mean Squared Residual')
ax2.scatter(damping_range, coeff_resid_list, label='Coefficient Mean Squared Residual')
ax1.set_xlabel('Damping Factor')
ax1.set_ylabel('Mean Squared Residual')
ax2.set_ylabel('Coefficient Mean Squared Residual')

ax1.legend(loc="upper left")
ax2.legend(loc="upper right")

# ax.set_title('Mean Squared Residual vs Damping Factor')
plt.savefig(f'{dest_folder}/RSSvsDamp.png')

# fig, ax = plt.subplots()

# ax.scatter(damping_range, r_squared_list)
# ax.set_xlabel('Damping Factor')
# ax.set_ylabel(r'$r^2$')
# ax.set_title(r'$r^2$ vs Damping Factor')
# plt.savefig(f'{dest_folder}/r2vsDamp.png')
