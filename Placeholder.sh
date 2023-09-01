


# *********************************************************************
# CENTRAL FILE OF COMMAND STRUCTURE FOR PROGRAM
#
# *********************************************************************

source Startup.sh
source Run_CUBIT.sh
yes | cp -rf geometry.exo geometry_base.exo
##Get the base version and run Pylith
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
yes | cp -rf UpperSurf_t0000000.vtk UpperSurf_base.vtk
cd ..


##Now modify thickness run Pylith to get Strain
python Modify_Thickness.py
yes | cp -rf Thickness_Profile.txt $cwd/Results/Thickness_Profile_true.txt
source Run_Pylith_dummy.sh
source shell/fix_path.sh
cd output 
yes | cp -rf All_Nodes_t0000000.vtk $cwd
yes | cp -rf All_Nodes_t0000000.vtk $cwd/Results/All_Nodes_true.vtk
yes | cp -rf LeftSurf_t0000000.vtk $cwd
yes | cp -rf RightSurf_t0000000.vtk $cwd
cd $cwd
python Make_ForceFile.py
yes | cp -rf forces.spatialdb spatialdb/forces.spatialdb
source Run_Pylith.sh
source shell/fix_path.sh
cd output
yes | cp -rf UpperSurf_t0000000.vtk $cwd
yes | cp -rf UpperSurf_t0000000.vtk UpperSurf_true.vtk
yes | cp -rf UpperSurf_t0000000.vtk $cwd/Results/UpperSurf_true.vtk
##Get Strain Profile and Create new geometry
yes | cp -rf UpperSurf_t0000000.vtk UpperSurf_thick.vtk
python Process_Thick.py
yes | cp -rf Strain_Profile.txt $cwd
yes | cp -rf Strain_Profile.txt $cwd/Results/Strain_Profile_true.txt
cd $cwd
yes | cp -rf geometry_base.exo geometry.exo
python Modify_Thickness_Strain.py







