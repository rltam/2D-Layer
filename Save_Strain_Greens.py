import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
from utility import Parse_vtk_1D
import math
import shutil

# finame="UpperSurf_t0000000"
# NodesData,DispData=Parse_vtk_1D.main(finame)[0],Parse_vtk_1D.main(finame)[2]

# DispData = DispData[DispData[:,0].argsort()]
# NodesData = NodesData[NodesData[:,0].argsort()]


# Displacementx = DispData[:,0]
# Displacementy = DispData[:,1]
# Locx = NodesData[:,0]
# Disp_Results_base=np.zeros((len(DispData[:,0]),3))
# Disp_Results_base[:,0]=Locx 
# Disp_Results_base[:,1]=Displacementx 
# Disp_Results_base[:,2]=Displacementy


# strain=np.zeros(len(Locx))
# for i in range(len(Locx)-1):
# 	strain[i]=(Disp_Results_base[i+1,1]-Disp_Results_base[i,1])/(Disp_Results_base[i+1,0]-Disp_Results_base[i,0])


try:
	iter_num = sys.argv[1]
except IndexError:
	iter_num = 40


if int(iter_num) < 0:
    shutil.copy('Diff_Strain_Profile.txt', 'GreenFunc/x_strain_unperturbed.txt')

else:
    shutil.copy('Diff_Strain_Profile.txt', f'GreenFunc/x_strain_{iter_num}.txt')
