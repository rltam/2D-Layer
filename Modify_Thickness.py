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
import scipy
from scipy import interpolate
from scipy.interpolate import BSpline, splrep
from Basis_Funcs import Spline_Basis

try:
	iter_num = int(sys.argv[4])
except IndexError:
	iter_num = 0

try:
	n_spline = int(sys.argv[1])
except IndexError:
	n_spline = 25

try:
	deg = int(sys.argv[2])
except IndexError:
	deg = 3

try:
	perc_amp = float(sys.argv[3])/100
except IndexError:
	perc_amp = 0.01

try:
	osc = int(sys.argv[5]) * 5
except IndexError:
	osc = 0

Edge_Height=np.loadtxt("Edge_Height.txt")
Edge_Length=np.loadtxt("Edge_Length.txt")


# ----------------------------------------------------------------------
# Get coordinates of points from ExodusII file.
exodus = netCDF4.Dataset(finame, 'a')
coords = np.empty((len(exodus.variables['coordx'][:]),2))
coords[:,0] = exodus.variables['coordx'][:]
coords[:,1] = exodus.variables['coordy'][:]

# ----------------------------------------------------------------------

coordsold=coords

##

basis = Spline_Basis(n_spline, deg)

# def Get_fun(coordsx,coordsy,maxcoordsx,maxcoordsy):
# 	'''
# 	Used to generate Green's matrix from basis
# 	'''
# 	amp =Edge_Height * perc_amp
# 	coeffs = np.zeros(len(basis))
# 	coeffs[int(iter_num)] = 1
# 	np.savetxt('Coefficients.txt', coeffs)

# 	if int(iter_num) < 0:
# 		Val = 0
# 	elif int(iter_num) == len(basis)-1:
# 		func = basis[int(iter_num)]
# 		if coordsx >= 49999:
# 			Val = 1 * amp
# 		else:
# 			Val = max(0, func(coordsx)) * amp
# 	else:
# 		func = basis[int(iter_num)]
# 		Val = max(0, func(coordsx)) * amp

# 	return Val

# def Get_fun(coordsx,coordsy,maxcoordsx,maxcoordsy):
# 	'''
# 	sin func
# 	'''
# 	# amp = Edge_Height * perc_amp
# 	osc = 7
# 	Wavelength = maxcoordsx *2/osc
# 	omega=2*np.pi/(Wavelength)
# 	Val= 20*np.sin(omega*(coordsx))
# 	file = open('Coefficients.txt', 'w')
# 	file.write(f'Sine Function with freq = {omega}')
# 	file.close()
# 	return Val

# def Get_fun(coordsx,coordsy,maxcoordsx,maxcoordsy):
# 	'''
# 	ramp
# 	'''
# 	amp = 0.001
# 	Val = amp * np.abs(coordsx)
# 	maxheight = amp * maxcoordsx
# 	file = open('Coefficients.txt', 'w')
# 	file.write(f'Linear Ramp with max height = {maxheight}')
# 	file.close()
# 	return Val

# def Get_fun(coordsx,coordsy,maxcoordsx,maxcoordsy):
# 	'''
# 	DC offset
# 	'''
# 	Val = Edge_Height * 0.01
# 	return Val

# def Get_fun(coordsx,coordsy,maxcoordsx,maxcoordsy):
# 	'''
# 	raised cosine
# 	'''
# 	wavelength = maxcoordsx
# 	omega = 2*np.pi/wavelength
# 	amp = perc_amp * Edge_Height / 2
# 	if np.abs(coordsx) <= wavelength/2:
# 		Val = amp * (np.cos(omega*coordsx)+1) * 5
# 	else:
# 		Val = 0
# 	return Val

# coeffs = np.zeros(n_spline)
# coeffs[7] = 0
# coeffs[10] = 1
# coeffs[17] = 10
# coeffs = np.loadtxt('Coefficients.txt')

coeffs = np.random.uniform(-10, 10, n_spline)

if osc > 0:
	coeffs = np.loadtxt('DerivedCoefficients.txt')

def Get_fun(coordsx,coordsy,maxcoordsx,maxcoordsy):
	'''
	superposition of basis fxns
	'''
	Val = 0
	amp = Edge_Height * perc_amp
	np.savetxt('Coefficients.txt', coeffs)

	for i in range(len(basis)):
		coeff = coeffs[i]
		func = basis[i]
		if i == len(basis)-1:
			if coordsx > 49999:
			# prevent last spline from dropping steeply to 0 at the end
 				Val += amp * coeff
			else:
				Val += max(0, func(coordsx)) * amp * coeff
		else:
			Val += max(0, func(coordsx)) * amp * coeff
	
	return Val+500



for i in range(len(coords[:,0])):
	if coords[i,1]< -Edge_Height/2+2:
		#print(1)
		coords[i,1]=coords[i,1]-Get_fun(coords[i,0],coords[i,1],np.amax(coords[:,0]),np.amax(coords[:,1]))

	elif coords[i,1]> Edge_Height/2-2:

		coords[i,1]=coords[i,1]+Get_fun(coords[i,0],coords[i,1],np.amax(coords[:,0]),np.amax(coords[:,1]))


	elif coords[i,1]> -Edge_Height/2+2 and coords[i,1]< Edge_Height/2-2:
		Old_Bot= -Edge_Height/2
		Old_Top= Edge_Height/2
		Loc_norm= (coords[i,1]-Old_Bot)/(Old_Top-Old_Bot)
		New_Bot= Old_Bot-Get_fun(coords[i,0],Old_Bot,np.amax(coords[:,0]),np.amax(coords[:,1]))
		New_Top= Old_Top+Get_fun(coords[i,0],Old_Bot,np.amax(coords[:,0]),np.amax(coords[:,1]))
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
	Thickness[i]=Edge_Height+2*Get_fun(coords[i,0],coords[i,1],np.amax(coords[:,0]),np.amax(coords[:,1]))#+np.abs(coords[i,0])/np.amax(coords[:,0])*300#*0

Profile=np.zeros((len(Thickness),2))
Profile[:,0]=coords[:,0]
Profile[:,1]=Thickness


np.savetxt("Thickness_Profile.txt",Profile)

exodus.close()
print('\nNew Thickness Saved!\n')
