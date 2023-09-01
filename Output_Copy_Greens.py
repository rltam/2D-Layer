import sys
import shutil
import os
import datetime
from datetime import datetime

txtfiles = ['Displacement_Profile', 'Strain_Profile', 'Thickness_Profile', 'Diff_Strain_Profile', 'Coefficients']
vtkfiles = ['UpperSurf', 'LeftSurf', 'RightSurf', 'All_Nodes']

try:
	iter_num = int(sys.argv[1])
except IndexError:
	iter_num = 0

try:
	n_spline = int(sys.argv[2])
except IndexError:
	n_spline = 25

try:
	deg = int(sys.argv[3])
except IndexError:
	deg = 3

try:
	amp = int(sys.argv[4])
except IndexError:
	amp = 1

try:
    iteration = int(sys.argv[5])
except IndexError:
    iteration = 0

try:
    folder = sys.argv[6]
except IndexError:
    folder = 'misc'


if iteration > 0:
    try:
        os.mkdir(f'/home/rltam/Documents/RunOutputFiles/{folder}/{iteration}_Greens({n_spline},{deg})x{amp}')
        big_folder = f'/home/rltam/Documents/RunOutputFiles/{folder}/{iteration}_Greens({n_spline},{deg})x{amp}'
    except FileExistsError:
        big_folder = f'/home/rltam/Documents/RunOutputFiles/{folder}/{iteration}_Greens({n_spline},{deg})x{amp}'
    
    os.mkdir(f'/home/rltam/Documents/RunOutputFiles/{folder}/{iteration}_Greens({n_spline},{deg})x{amp}/Spline{iter_num}')
    dest_folder = f'{big_folder}/Spline{iter_num}'

else:
    try:
        os.mkdir(f'/home/rltam/Documents/RunOutputFiles/Greens({n_spline},{deg})x{amp}')
        big_folder = f'/home/rltam/Documents/RunOutputFiles/Greens({n_spline},{deg})x{amp}'
    except FileExistsError:
        big_folder = f'/home/rltam/Documents/RunOutputFiles/Greens({n_spline},{deg})x{amp}'

    os.mkdir(f'/home/rltam/Documents/RunOutputFiles/Greens({n_spline},{deg})x{amp}/Spline{iter_num}')
    dest_folder = f'{big_folder}/Spline{iter_num}'


for file in txtfiles:
    txtSplit = file.split('_')
    word = txtSplit[0]
    shutil.copy(f'{file}.txt', f'{dest_folder}/{file}.txt')

for file in vtkfiles:
    txtSplit = file.split('_')
    word = txtSplit[0]
    shutil.copy(f'{file}_t0000000.vtk', f'{dest_folder}/{file}.vtk')

shutil.copy('ThicknessAndStrain.png', f'{dest_folder}/ThicknessAndStrain.png')
shutil.copy('ThicknessAndStrainDiff.png', f'{dest_folder}/ThicknessAndStrainDiff.png')
