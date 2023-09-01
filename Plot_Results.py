##Plotting Results....

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
from utility import Parse_vtk_1D

try:
	iter_num = sys.argv[1]
except IndexError:
	iter_num = 0

finame="UpperSurf_t0000000"
NodesData,DispData=Parse_vtk_1D.main(finame)[0],Parse_vtk_1D.main(finame)[2]


Displacementx = DispData[:,0]
Displacementy = DispData[:,1]
Locx = NodesData[:,0]

Disp_Results=np.zeros((len(DispData[:,0]),3))
Disp_Results[:,0]=Locx 
Disp_Results[:,1]=Displacementx 
Disp_Results[:,2]=Displacementy
np.savetxt('Displacement_Profile.txt',Disp_Results)

Thickness=np.loadtxt('Thickness_Profile.txt')

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()  # Create a twinned axis

ax1.scatter(Thickness[:, 0], Thickness[:, 1], c='r', s=10, label="Thickness")
ax2.scatter(Disp_Results[:, 0], Disp_Results[:, 1], c='b', s=10, label="Displacement X")
ax2.scatter(Disp_Results[:, 0], Disp_Results[:, 2], c='g', s=10, label="Displacement Y")

ax1.set_xlabel("Distance Along Block (m)")
ax1.set_ylabel("Thickness (m)")
ax2.set_ylabel("Displacement (m)")

ax1.legend(loc="upper left")
ax2.legend(loc="upper right")

plt.show()
# plt.savefig(f'/home/rltam/2D_Layer-main/plot_func{iter_num}.png')








