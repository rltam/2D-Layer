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
txtfiles = ['Displacement_Profile', 'Strain_Profile', 'Thickness_Profile', 'Diff_Strain_Profile']
vtkfiles = ['UpperSurf', 'LeftSurf', 'RightSurf', 'All_Nodes']

now = datetime.now()
dt_string = now.strftime('%m.%d.%Y_%H:%M:%S')

os.mkdir(f'/home/rltam/Documents/RunOutputFiles/{dt_string}')
os.mkdir(f'/home/rltam/Documents/RunOutputFiles/{dt_string}/RunData')
os.mkdir(f'/home/rltam/Documents/RunOutputFiles/{dt_string}/InverseOutputs')


dest_folder = f'/home/rltam/Documents/RunOutputFiles/{dt_string}'

for file in txtfiles:
    txtSplit = file.split('_')
    word = txtSplit[0]
    shutil.copy(f'output/{file}.txt', f'{dest_folder}/RunData/{word}_{dt_string}.txt')

for file in vtkfiles:
    txtSplit = file.split('_')
    word = txtSplit[0]
    shutil.copy(f'output/{file}_t0000000.vtk', f'{dest_folder}/RunData/{word}_{dt_string}.vtk')

shutil.copy('Coefficients.txt', f'{dest_folder}/RunData/Coefficients_{dt_string}.txt')

# everything else

# strain = np.loadtxt('output/Strain_Profile.txt')
# base = np.loadtxt('GreenFunc/x_strain_unperturbed.txt')
Edge_Height=np.loadtxt("Edge_Height.txt")
Edge_Length=np.loadtxt("Edge_Length.txt")
true_thick = np.loadtxt('output/Thickness_Profile.txt')
input_coeffs = np.loadtxt('Coefficients.txt')

frac_strain = np.loadtxt('Diff_Strain_Profile.txt')

G30 = np.loadtxt('GreenFunc/matrix30.txt')
G30_offset = np.loadtxt('GreenFunc/matrixOffset30.txt')
G50 = np.loadtxt('GreenFunc/matrix50.txt')
G50_offset = np.loadtxt('GreenFunc/matrixOffset50.txt')


## cut out edges:
G30 = G30[:-1]
G30_offset = G30_offset[:-1]
G50 = G50[:-1]
G50_offset = G50_offset[:-1]

frac_strain_cut = frac_strain[:-1]

basis30 = Spline_Basis(30, 3)
basis50 = Spline_Basis(50, 3)

true_thick = true_thick[true_thick[:,0].argsort()]
x = true_thick[:,0]
true_thickness = [true_thick[:,1][i] - 1000 for i in range(len(true_thick[:,1]))]

# damping_range = np.linspace(5e-8, 5e-7, 10)
damping_range = np.logspace(-9, -6, 7)

resid_list_30 = []
resid_list_30_offset = []
resid_list_50 = []
resid_list_50_offset = []


for damping in damping_range:
    print(f'Now using damping factor {damping}')

    ## do inverse
    inverse30 = scipy.sparse.linalg.lsqr(G30, frac_strain_cut, damp = damping, show=False)
    coeffs30 = inverse30[0]
    np.savetxt(f'{dest_folder}/InverseOutputs/InverseCoeffs30_{damping}', coeffs30)

    inverse30_offset = scipy.sparse.linalg.lsqr(G30_offset, frac_strain_cut, damp = damping, show=False)
    coeffs30_offset = inverse30_offset[0]
    np.savetxt(f'{dest_folder}/InverseOutputs/InverseCoeffs30Offset_{damping}', coeffs30_offset)

    inverse50 = scipy.sparse.linalg.lsqr(G50, frac_strain_cut, damp = damping, show=False)
    coeffs50 = inverse50[0]
    np.savetxt(f'{dest_folder}/InverseOutputs/InverseCoeffs50_{damping}', coeffs50)

    inverse50_offset = scipy.sparse.linalg.lsqr(G50_offset, frac_strain_cut, damp = damping, show=False)
    coeffs50_offset = inverse50_offset[0]
    np.savetxt(f'{dest_folder}/InverseOutputs/InverseCoeffs50Offset_{damping}', coeffs50_offset)
    
    ## derive thickness profiles
    derived_prof_30 = []
    derived_prof_30_offset = []
    derived_prof_50 = []
    derived_prof_50_offset = []

    for i in range(len(x)):
        val30 = 0
        val30_offset = 0
        val50 = 0
        val50_offset = 0

        for j in range(len(basis30)):
            base = basis30[j]
            if j == len(basis30) - 1:
                if x[i] >= 49999:
                    val30 += 1 * Edge_Height * 0.01 * coeffs30[j]
                    val30_offset += 1 * Edge_Height * 0.01 * coeffs30_offset[j]
                else:
                    val30 += Edge_Height * 0.01*coeffs30[j]*max(base(x[i]), 0)
                    val30_offset += Edge_Height * 0.01 * coeffs30_offset[j] * max(base(x[i]), 0)
            else:
                val30 += Edge_Height * 0.01 * coeffs30[j] * max(base(x[i]), 0)
                val30_offset += Edge_Height * 0.01 * coeffs30_offset[j] * max(base(x[i]), 0)
        
        val30_offset += Edge_Height * 0.01 * coeffs30_offset[-1]   
        
        derived_prof_30.append(val30)
        derived_prof_30_offset.append(val30_offset)


        for j in range(len(basis50)):
            base = basis50[j]
            if j == len(basis50) - 1:
                if x[i] >= Edge_Length/2 - 1:
                    val50 += 1 * Edge_Height * 0.01 * coeffs50[j]
                    val50_offset += 1 * Edge_Height * 0.01 * coeffs50_offset[j]
                else:
                    val50 += Edge_Height * 0.01*coeffs50[j]*max(base(x[i]), 0)
                    val50_offset += Edge_Height * 0.01 * coeffs50_offset[j] * max(base(x[i]), 0)
            else:
                val50 += Edge_Height * 0.01 * coeffs50[j] * max(base(x[i]), 0)
                val50_offset += Edge_Height * 0.01 * coeffs50_offset[j] * max(base(x[i]), 0)
        
        val50_offset += Edge_Height * 0.01 * coeffs50_offset[-1]   
        
        derived_prof_50.append(val50)
        derived_prof_50_offset.append(val50_offset)

    ## residuals 

    diffs30 = np.zeros(len(x))
    diffs30_offset = np.zeros(len(x))
    diffs50 = np.zeros(len(x))
    diffs50_offset = np.zeros(len(x))

    for i in range(len(x)):
        diffs30[i] = derived_prof_30[i] - true_thickness[i]
        diffs30_offset[i] = derived_prof_30_offset[i] - true_thickness[i]
        diffs50[i] = derived_prof_50[i] - true_thickness[i]
        diffs50_offset[i] = derived_prof_50_offset[i] - true_thickness[i]

    diffs_squared_30 = [(diffs30[i])**2 for i in range(len(x))]
    diffs_squared_30_offset = [(diffs30_offset[i])**2 for i in range(len(x))]
    diffs_squared_50 = [(diffs50[i])**2 for i in range(len(x))]
    diffs_squared_50_offset = [(diffs50_offset[i])**2 for i in range(len(x))]

    meanRS30 = np.mean(diffs_squared_30)
    meanRS30_offset = np.mean(diffs_squared_30_offset)
    meanRS50 = np.mean(diffs_squared_50)
    meanRS50_offset = np.mean(diffs_squared_50_offset)

    all_diffs = np.zeros((len(diffs30), 4))
    all_diffs[:,0] = diffs30
    all_diffs[:,1] = diffs30_offset
    all_diffs[:,2] = diffs50
    all_diffs[:,3] = diffs50_offset
    np.savetxt(f'{dest_folder}/InverseOutputs/Residuals_{damping}.txt',all_diffs)

    all_diffs_squared = np.zeros((len(diffs_squared_30), 4))
    all_diffs_squared[:,0] = diffs_squared_30
    all_diffs_squared[:,1] = diffs_squared_30_offset
    all_diffs_squared[:,2] = diffs_squared_50
    all_diffs_squared[:,3] = diffs_squared_50_offset

    np.savetxt(f'{dest_folder}/InverseOutputs/Residuals2_{damping}.txt',all_diffs_squared)
    
    resid_list_30.append(meanRS30)
    resid_list_30_offset.append(meanRS30_offset)
    resid_list_50.append(meanRS50)
    resid_list_50_offset.append(meanRS50_offset)


    # make plot
    fig, ax = plt.subplots(3,2, figsize=(20, 22))

    ax[0][0].plot(x, derived_prof_30, c = 'tomato', label = 'Derived Profile, 30 Basis Functions')
    ax[0][0].plot(x, derived_prof_30_offset, c = 'darkorange', label = 'Derived Profile, 30 Basis Functions + Offset')
    ax[0][0].plot(x, derived_prof_50, c = 'lightseagreen', label = 'Derived Profile, 50 Basis Functions')
    ax[0][0].plot(x, derived_prof_50_offset, c = 'royalblue', label = 'Derived Profile, 50 Basis Functions + Offset')
    ax[0][0].plot(true_thick[:,0], [true_thick[:,1][i]-1000 for i in range(len(true_thick[:,1]))], c = 'k', linewidth=2, label = 'True Profile')

    ax[0][0].legend()
    ax[0][0].set_xlabel("Distance Along Block (m)")
    ax[0][0].set_ylabel("Thickness Variation (m)")
    ax[0][0].set_title('True and Derived Thickness Profiles')

    ax[0][1].scatter(x, diffs30, s = 5, c = 'tomato', label='30 Basis Functions')
    ax[0][1].scatter(x, diffs30_offset, s = 5, c = 'darkorange', label='30 Basis Functions + Offset')
    ax[0][1].scatter(x, diffs50, s = 5, c = 'lightseagreen', label='50 Basis Functions')
    ax[0][1].scatter(x, diffs50_offset, s = 5, c = 'royalblue', label='50 Basis Functions + Offset')
    ax[0][1].grid()
    ax[0][1].set_xlabel("Distance Along Block (m)")
    ax[0][1].set_ylabel('Residual')
    ax[0][1].set_title('Residuals')
    ax[0][1].legend()

    hist(diffs_squared_30, bins = 30, ax = ax[1][0], color = 'tomato', label=f'Mean Squared Residual = {meanRS30}')
    hist(diffs_squared_30_offset, bins = 30, ax = ax[1][1], color = 'darkorange',  label=f'Mean Squared Residual = {meanRS30_offset}')
    hist(diffs_squared_50, bins = 30, ax = ax[2][0], color = 'lightseagreen',  label=f'Mean Squared Residual = {meanRS50}')
    hist(diffs_squared_50_offset, bins = 30, ax = ax[2][1], color = 'royalblue',  label=f'Mean Squared Residual = {meanRS50_offset}')


    ax[1][0].set_xlabel('Squared Residual')
    ax[1][0].set_ylabel('Number of Data Points') 
    ax[1][0].set_title('30 Basis Functions')
    ax[1][0].legend()

    ax[1][1].set_xlabel('Squared Residual')
    ax[1][1].set_ylabel('Number of Data Points') 
    ax[1][1].set_title('30 Basis Functions + Offset')
    ax[1][1].legend()

    ax[2][0].set_xlabel('Squared Residual')
    ax[2][0].set_ylabel('Number of Data Points') 
    ax[2][0].set_title('50 Basis Functions')
    ax[2][0].legend()

    ax[2][1].set_xlabel('Squared Residual')
    ax[2][1].set_ylabel('Number of Data Points') 
    ax[2][1].set_title('50 Basis Functions + Offset')
    ax[2][1].legend()

    fig.suptitle(f'Profile and Residual Histogram for Inverse Method, damp = {damping}')

    # plt.show()
    plt.savefig(f'{dest_folder}/plots_{damping}.png')


mean_data = np.zeros((len(damping_range), 5))

mean_data[:,0] = damping_range
mean_data[:,1] = resid_list_30
mean_data[:,2] = resid_list_30_offset
mean_data[:,3] = resid_list_50
mean_data[:,4] = resid_list_50_offset

np.savetxt(f'{dest_folder}/InverseOutputs/MeanSquaredResiduals.txt', mean_data)

fig, ax = plt.subplots(figsize=(12, 8))

ax.scatter(damping_range, resid_list_30, c='tomato', label='30 Basis Functions')
ax.scatter(damping_range, resid_list_30_offset, c='darkorange', label='30 Basis Functions + Offset')
ax.scatter(damping_range, resid_list_30, c='lightseagreen', label='50 Basis Functions')
ax.scatter(damping_range, resid_list_30_offset, c='royalblue', label='50 Basis Functions + Offset')

ax.set_xlabel('Damping Factor')
ax.set_ylabel('Mean Squared Residual')
ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid()
ax.set_title('Mean Squared Residual vs Damping Factor')
plt.savefig(f'{dest_folder}/RSSvsDamp.png')
