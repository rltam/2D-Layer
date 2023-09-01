


# *********************************************************************
# CENTRAL FILE OF COMMAND STRUCTURE FOR PROGRAM
#
# *********************************************************************


source Startup.sh
source Run_CUBIT.sh
python Modify_Thickness.py $i
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
python Plot_Results.py $i
python Save_Strain.py $i

yes | cp -rf Displacement_Profile.txt $cwd/output
yes | cp -rf Thickness_Profile.txt $cwd/output
yes | cp -rf Strain_Profile.txt $cwd/output
yes | cp -rf Diff_Strain_Profile.txt $cwd/output

# python Output_Copy.py

