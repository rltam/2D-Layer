


# *********************************************************************
# CENTRAL FILE OF COMMAND STRUCTURE FOR PROGRAM
#
# *********************************************************************

#From Greens: $1 = number of splines , $2 = degree , $3 = amplitude , $4 = perturbation iteration #, $5 = solver iteration number, $6 = folder

source Startup.sh
source Run_CUBIT.sh
python Modify_Thickness_Greens.py $1 $2 $3 $4
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
# python Plot_Results.py $i
# python Save_Strain.py $i


# yes | cp -rf Displacement_Profile.txt $cwd/output
# yes | cp -rf Thickness_Profile.txt $cwd/output
# yes | cp -rf Strain_Profile.txt $cwd/output
# yes | cp -rf Diff_Strain_Profile.txt $cwd/output

python Output_Copy_Greens.py $4 $1 $2 $3 $5 $6

