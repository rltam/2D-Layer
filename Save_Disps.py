import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
from utility import Parse_vtk_1D
import math


iter_num = int(np.loadtxt('Iteration.txt'))

disp_file = open('output/Displacement_Profile.txt', 'r')
lines = disp_file.readlines()
disp_file.close()

x_block = []
x_disp = []
y_disp = []

for line in lines:
    line.strip()
    nums = line.split()
    x_block.append(float(nums[0]))
    x_disp.append(float(nums[1]))
    y_disp.append(float(nums[2]))

if iter_num < 0:
    np.savetxt('GreenFunc/x_disp_unperturbed.txt', x_disp)
else:
    unperturbed = np.loadtxt('GreenFunc/x_disp_unperturbed.txt')
    frac_disp = []
    for j in range(len(unperturbed)):
        frac = x_disp[j]/unperturbed[j]
        frac_disp.append(frac)
    np.savetxt(f'GreenFunc/x_disp_{iter_num}.txt', frac_disp)

