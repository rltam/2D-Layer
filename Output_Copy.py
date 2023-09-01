import sys
import shutil
import os
import datetime
from datetime import datetime

txtfiles = ['Displacement_Profile', 'Strain_Profile', 'Thickness_Profile', 'Diff_Strain_Profile', 'Coefficients']
vtkfiles = ['UpperSurf', 'LeftSurf', 'RightSurf', 'All_Nodes']

folder = sys.argv[1]

try:
    os.mkdir(f'/home/rltam/Documents/RunOutputFiles/{folder}')
except FileExistsError:
    print('folder exists')

try:
    os.mkdir(f'/home/rltam/Documents/RunOutputFiles/{folder}/True_Data')
    dest_folder = f'/home/rltam/Documents/RunOutputFiles/{folder}/True_Data'
except FileExistsError:
    dest_folder = f'/home/rltam/Documents/RunOutputFiles/{folder}/True_Data'

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
