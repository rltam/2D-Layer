##Plotting Results....

import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
from utility import Parse_vtk_1D
import shutil

# try:
#     folder = sys.argv[1]
# except IndexError:
#     folder = 'warre'

# try:
#     os.mkdir(f'/home/rltam/Documents/RunOutputFiles/{folder}')
#     dest_folder = f'/home/rltam/Documents/RunOutputFiles/{folder}'
# except FileExistsError:
#     dest_folder = f'/home/rltam/Documents/RunOutputFiles/{folder}'
    
# finame="UpperSurf_t0000000"
# NodesData,DispData=Parse_vtk_1D.main(finame)[0],Parse_vtk_1D.main(finame)[2]


# Displacementx = DispData[:,0]
# Displacementy = DispData[:,1]
# Locx = NodesData[:,0]

# Disp_Results=np.zeros((len(DispData[:,0]),3))
# Disp_Results[:,0]=Locx 
# Disp_Results[:,1]=Displacementx 
# Disp_Results[:,2]=Displacementy
# np.savetxt('Displacement_Profile.txt',Disp_Results)

Thickness=np.loadtxt('Thickness_Profile.txt')
Thickness = Thickness[Thickness[:,0].argsort()]

Strain_Diff = np.loadtxt('Diff_Strain_Profile.txt')[:-1]
Strain = np.loadtxt('Strain_Profile.txt')[:-1]

fig, ax1 = plt.subplots(figsize = (10,7))
ax2 = ax1.twinx()  # Create a twinned axis
# ax1.scatter(Thickness[:, 0], Thickness[:, 1], c='r', s=5, label="Thickness")
# ax2.scatter(Strain_Diff[:,0], Strain_Diff[:,1], c='b', s=5, label = "Strain Difference")
ax1.plot(Thickness[:, 0], Thickness[:,1], c='r', label="Thickness")
ax2.plot(Strain_Diff[:,0], Strain_Diff[:,1], c='b', label = "Strain Difference")
ax1.set_xlabel("Distance Along Block (m)")
ax1.set_ylabel("Thickness (m)")
ax2.set_ylabel("Strain Diffence")
ax1.legend(loc="upper left")
ax2.legend(loc="upper right")
# ax2.set_ylim(-1.1e-6, -1e-6)
ax1.set_title('Thickness Profile and Strain Difference from Uniform Thickness Case')
plt.savefig('ThicknessAndStrainDiff.png')
# shutil.copy('ThicknessAndStrainDiff.png', f'{dest_folder}/ThicknessAndStrainDiff.png')
plt.close()

fig, ax1 = plt.subplots(figsize = (10,7))
ax2 = ax1.twinx()  # Create a twinned axis
# ax1.scatter(Thickness[:, 0], Thickness[:, 1], c='r', s=5, label="Thickness")
# ax2.scatter(Strain[:,0], Strain[:,1], c='b', s=5, label = "Strain")
ax1.plot(Thickness[:, 0], Thickness[:, 1], c='r', label="Thickness")
ax2.plot(Strain[:,0], Strain[:,1], c='b', label = "Strain")
ax1.set_xlabel("Distance Along Block (m)")
ax1.set_ylabel("Thickness (m)")
ax2.set_ylabel("Strain")
# ax2.set_ylim(5.4e-6, 5.45e-6)
ax1.legend(loc="upper left")
ax2.legend(loc="upper right")

ax1.set_title('Thickness and Strain Profiles Along Block')

plt.savefig('ThicknessAndStrain.png')
# shutil.copy('ThicknessAndStrain.png', f'{dest_folder}/ThicknessAndStrain.png')
plt.close()

