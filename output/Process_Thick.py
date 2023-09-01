

##Plotting Results....

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
from utility import Parse_vtk_1D


finame="UpperSurf_base"
NodesData,DispData=Parse_vtk_1D.main(finame)[0],Parse_vtk_1D.main(finame)[2]


Displacementx = DispData[:,0]
Displacementy = DispData[:,1]
Locx = NodesData[:,0]
Disp_Results_base=np.zeros((len(DispData[:,0]),3))
Disp_Results_base[:,0]=Locx 
Disp_Results_base[:,1]=Displacementx 
Disp_Results_base[:,2]=Displacementy


strainxx_base=np.zeros(len(Locx))
for i in range(len(Locx)-1):
	strainxx_base[i]=(Disp_Results_base[i+1,1]-Disp_Results_base[i,1])/(Disp_Results_base[i+1,0]-Disp_Results_base[i,0])


finame="UpperSurf_thick"
NodesData,DispData=Parse_vtk_1D.main(finame)[0],Parse_vtk_1D.main(finame)[2]


Displacementx = DispData[:,0]
Displacementy = DispData[:,1]
Locx = NodesData[:,0]
Disp_Results=np.zeros((len(DispData[:,0]),3))
Disp_Results[:,0]=Locx 
Disp_Results[:,1]=Displacementx 
Disp_Results[:,2]=Displacementy


strainxx=np.zeros(len(Locx))
for i in range(len(Locx)-1):
	strainxx[i]=(Disp_Results[i+1,1]-Disp_Results[i,1])/(Disp_Results[i+1,0]-Disp_Results[i,0])



Strain_ratio=strainxx/strainxx_base

#print(Strain_ratio)
#plt.scatter(Locx,Strain_ratio)
#plt.show()

Strain_Profile=np.zeros((len(Strain_ratio),2))
Strain_Profile[:,0]= Locx 
Strain_Profile[:,1]= np.nan_to_num(Strain_ratio,nan=1.0)

np.savetxt('Strain_Profile.txt',Strain_Profile)







