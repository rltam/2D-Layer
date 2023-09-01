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
import shutil

# np.savetxt('Iteration.txt',[int(-1)])

# # subprocess.Popen("ls", stdout=subprocess.PIPE, shell=True)

# subprocess.run(["/home/rltam/2DLayer/Go.sh"])

# disp_file_unperturbed = open('output/Displacement_Profile.txt', 'r')
# lines = disp_file_unperturbed.readlines()
# disp_file_unperturbed.close()

# x_block = []
# x_disp_unperturbed = []

# for line in lines:
#     line.strip()
#     nums = line.split()
#     x_block.append(float(nums[0]))
#     x_disp_unperturbed.append(float(nums[1]))

try:
	n_spline = int(sys.argv[1])
except IndexError:
	n_spline = 25

try:
	deg = int(sys.argv[2])
except IndexError:
	deg = 3


try:
    amp = int(sys.argv[3])
except IndexError:
    amp = 1

try:
    iteration = int(sys.argv[4])
except IndexError:
    iteration = -1


try:
    folder = sys.argv[5]
except IndexError:
    folder = 'misc'

basis = Spline_Basis(n_spline, deg)

# for i in range(len(basis)):
#     np.savetxt('Iteration.txt',[i])

#     # os.system("sh Go.sh")

#     disp_file_perturbed = open('output/Displacement_Profile.txt', 'r')
#     lines = disp_file_perturbed.readlines()
#     disp_file_perturbed.close()
#     x_disp_perturbed = []
#     for line in lines:
#         line.strip()
#         nums = line.split()
#         x_disp_perturbed.append(float(nums[1]))

#     frac_disp = []
#     for j in range(len(x_disp_perturbed)):
#         frac = x_disp_perturbed[i]/x_disp_unperturbed[i]
#         frac_disp.append(frac)

#     np.savetxt(f'GreenFunc/x_disp_{i}.txt', frac_disp)

## make a file with the green's function

# file = open(f'GreenFunc/matrix({n_spline},{deg}).txt', 'w')
# file.close()

dest_folder = f'/home/rltam/Documents/RunOutputFiles/{folder}/{iteration}_Greens({n_spline},{deg})x{amp}'

unperturbed = np.loadtxt('GreenFunc/x_strain_unperturbed.txt')

G = np.zeros((len(unperturbed),len(basis)))

for i in range(len(basis)):
    col = np.loadtxt(f'GreenFunc/x_strain_{i}.txt')
    G[:,i] = col[:,1]

np.savetxt(f'GreenFunc/matrix({n_spline},{deg})x{amp}_{iteration}.txt', G)
shutil.copy(f'GreenFunc/matrix({n_spline},{deg})x{amp}_{iteration}.txt', f'{dest_folder}/matrix({n_spline},{deg})x{amp}_{iteration}.txt')