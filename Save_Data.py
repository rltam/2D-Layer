import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
from utility import Parse_vtk_1D
import math

try:
	iteration = int(sys.argv[1])
except IndexError:
	iteration = -1

finame="UpperSurf_t0000000"
NodesData,DispData=Parse_vtk_1D.main(finame)[0],Parse_vtk_1D.main(finame)[2]

DispData = DispData[DispData[:,0].argsort()]
NodesData = NodesData[NodesData[:,0].argsort()]


Displacementx = DispData[:,0]
Displacementy = DispData[:,1]
Locx = NodesData[:,0]
Disp_Results_base=np.zeros((len(DispData[:,0]),3))
Disp_Results_base[:,0]=Locx 
Disp_Results_base[:,1]=Displacementx 
Disp_Results_base[:,2]=Displacementy


np.savetxt('Displacement_Profile.txt', Disp_Results_base)

strain=np.zeros(len(Locx))
Location=np.zeros(len(Locx))
for i in range(len(Locx)-1):
	strain[i]=(Disp_Results_base[i+1,1]-Disp_Results_base[i,1])/(Disp_Results_base[i+1,0]-Disp_Results_base[i,0])

for i in range(len(Locx)-1):
	Location[i]=(Disp_Results_base[i+1,0]-Disp_Results_base[i,0])/2+Disp_Results_base[i,0]


save_strain = np.zeros((len(strain),2))
save_strain[:,0] = Locx
save_strain[:,1] = strain

np.savetxt('Strain_Profile.txt', save_strain)

if iteration > 0:
	unperturbed = np.loadtxt('GreenFunc/x_strain_unperturbed.txt')
elif iteration <= 0:
	unperturbed = np.loadtxt('GreenFunc/Uniform_Strain_Profile.txt')

diff_strain = [strain[i] - unperturbed[:,1][i] for i in range(len(unperturbed))]

# diff_strain = np.zeros(len(unperturbed))
# for i in range(len(diff_strain)):
# 	diff_strain[i] = strain[i] - unperturbed[i]


save_diff_strain = np.zeros((len(diff_strain),2))
# save_diff_strain[:,0] = Locx
# save_diff_strain[:,1] = diff_strain

for i in range(len(save_diff_strain)):
	save_diff_strain[i][0] = Locx[i]
	save_diff_strain[i][1] = diff_strain[i]

np.savetxt('Diff_Strain_Profile.txt', save_diff_strain)

# thickness = np.loadtxt('')







