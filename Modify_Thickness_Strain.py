##Steps: 1. Mesh a coarse-ish sphere in CUBIT. NO UNTANGLING OR SMOOTHING

## 2. Export the mesh and run it through this program

## 3. Send it back to cubit and do a finer mesh

## 4. Run Simulation!


finame = "geometry.exo"

import sys
import numpy as np
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
import pandas as pd
import pyshtools as pysh
from pyshtools import constants
from pyshtools import gravmag
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import sph_harm
import math



Edge_Height=np.loadtxt("Edge_Height.txt")

Strain_Profile=np.loadtxt('Strain_Profile.txt')

# ----------------------------------------------------------------------
# Get coordinates of points from ExodusII file.
exodus = netCDF4.Dataset(finame, 'a')
coords = np.empty((len(exodus.variables['coordx'][:]),2))
coords[:,0] = exodus.variables['coordx'][:]
coords[:,1] = exodus.variables['coordy'][:]

Strain_Profile_n=np.zeros((len(coords[:,0]),2))

##Need to sort

Strain_Profile=Strain_Profile[Strain_Profile[:, 0].argsort()]

Strain_Profile_n[:,0]=coords[:,0]
Strain_Profile_n=Strain_Profile_n[Strain_Profile_n[:, 0].argsort()]
Strain_Profile_n[:,1]=np.interp(Strain_Profile_n[:,0],Strain_Profile[:,0],Strain_Profile[:,1])



##Need to re-sort strain profile into original configuration

Newvec=np.zeros((len(Strain_Profile_n[:,0]),2))

for i in range(len(Strain_Profile_n[:,0])):
	Diff=np.abs(coords[i,0]-Strain_Profile_n[:,0])
	index=np.argmin(Diff)
	Newvec[i,0]=Strain_Profile_n[int(index),0]
	Newvec[i,1]=Strain_Profile_n[int(index),1]
	
Strain_Profile=Newvec


# ----------------------------------------------------------------------

coordsold=coords

##


def Get_fun(Strain_Profile_loc,Strain_Profile_val,coordsx,coordsy,maxcoordsx,maxcoordsy):



	Val=(Strain_Profile_val-1)*0.5*Edge_Height

	#print(Val)

	
	return Val	

	


for i in range(len(coords[:,0])):
	if coords[i,1]< -Edge_Height/2+2:
		#print(1)
		coords[i,1]=coords[i,1]-Get_fun(Strain_Profile[i,0],Strain_Profile[i,1],coords[i,0],coords[i,1],np.amax(coords[:,0]),np.amax(coords[:,1]))#-np.abs(coords[i,0])/np.amax(coords[:,0])*300#*0

	elif coords[i,1]> -Edge_Height/2+2:
		Old_Bot= -Edge_Height/2
		Old_Top= Edge_Height/2
		Loc_norm= (coords[i,1]-Old_Bot)/(Old_Top-Old_Bot)
		New_Bot= Old_Bot-Get_fun(Strain_Profile[i,0],Strain_Profile[i,1],coords[i,0],Old_Bot,np.amax(coords[:,0]),np.amax(coords[:,1]))
		New_Top= Edge_Height/2
		coords[i,1]= (New_Top-New_Bot)*Loc_norm+New_Bot

##
print(np.amax(np.abs(coords-coordsold)))
xcoordsnew=coords[:,0]
ycoordsnew=coords[:,1]

##Put it all back into the exodus file

exodus.variables['coordx'][:]=xcoordsnew
exodus.variables['coordy'][:]=ycoordsnew


##Save txt file of thickness....

Thickness=np.zeros(len(xcoordsnew))

for i in range(len(Thickness)):
	Thickness[i]=Edge_Height+Get_fun(Strain_Profile[i,0],Strain_Profile[i,1],coords[i,0],coords[i,1],np.amax(coords[:,0]),np.amax(coords[:,1]))#+np.abs(coords[i,0])/np.amax(coords[:,0])*300#*0
	


Profile=np.zeros((len(Thickness),2))
Profile[:,0]=coords[:,0]
Profile[:,1]=Thickness


np.savetxt("Thickness_Profile.txt",Profile)

exodus.close()
