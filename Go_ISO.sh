
# *********************************************************************
# CENTRAL FILE OF COMMAND STRUCTURE FOR PROGRAM
#
# *********************************************************************

# argument entry:
# 1. output folder name
# 2. number of splines in basis
# 3. degree
# 4. amplitude

j=10

source Startup.sh
source Run_CUBIT.sh
python Modify_Thickness.py $2 $3 $4 $j
source Run_Pylith_dummy.sh
source shell/fix_path.sh
cd output 
yes | cp -rf All_Nodes_t0000000.vtk $cwd
yes | cp -rf LeftSurf_t0000000.vtk $cwd
yes | cp -rf RightSurf_t0000000.vtk $cwd
cd $cwd
python Make_ForceFile.py
yes | cp -rf forces.spatialdb spatialdb/forces.spatialdb
source Run_Pylith.sh
source shell/fix_path.sh
cd output
yes | cp -rf UpperSurf_t0000000.vtk $cwd
cd ..
python Save_Data.py
python Plots_ISO.py

# yes | cp -rf Displacement_Profile.txt $cwd/output
# yes | cp -rf Thickness_Profile.txt $cwd/output
# yes | cp -rf Strain_Profile.txt $cwd/output
# yes | cp -rf Diff_Strain_Profile.txt $cwd/output

python Output_Copy.py $1

