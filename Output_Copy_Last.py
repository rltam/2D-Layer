import sys
import shutil
import os
import datetime
from datetime import datetime
import numpy as np

txtfiles = ['Displacement_Profile', 'Strain_Profile', 'Thickness_Profile', 'Diff_Strain_Profile', 'Coefficients']
vtkfiles = ['UpperSurf', 'LeftSurf', 'RightSurf', 'All_Nodes']

folder = sys.argv[1]

os.mkdir(f'/home/rltam/Documents/RunOutputFiles/{folder}/Last_Iteration')

dest_folder = f'/home/rltam/Documents/RunOutputFiles/{folder}/Last_Iteration'

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

obs_strain = np.loadtxt(f'/home/rltam/Documents/RunOutputFiles/{folder}/True_Data/Strain_Profile.txt')
rec_strain = np.loadtxt('Strain_Profile.txt')
diff_strain = np.zeros(len(obs_strain))
for i in range(len(diff_strain)):
    diff_strain[i] = obs_strain[i][1] - rec_strain[i][1]
np.savetxt(f'/home/rltam/Documents/RunOutputFiles/{folder}/Last_Iteration/Obs-Rec_Strain.txt', diff_strain)