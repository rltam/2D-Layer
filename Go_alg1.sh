


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
python Modify_Thickness_First.py
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
yes | cp -rf Thickness_Profile.txt $cwd/Results/Thickness_Profile_n0s.txt
yes | cp -rf Strain_Profile.txt Strain_Profile_Prev.txt



##Now do the first Run
source Run_Pylith_dummy.sh
source shell/fix_path.sh
cd output 
yes | cp -rf All_Nodes_t0000000.vtk $cwd
yes | cp -rf All_Nodes_t0000000.vtk $cwd/Results/All_Nodes_n0s.vtk
yes | cp -rf LeftSurf_t0000000.vtk $cwd
yes | cp -rf RightSurf_t0000000.vtk $cwd
cd $cwd
python Make_ForceFile.py
yes | cp -rf forces.spatialdb spatialdb/forces.spatialdb
source Run_Pylith.sh
source shell/fix_path.sh
cd output
yes | cp -rf UpperSurf_t0000000.vtk $cwd
yes | cp -rf UpperSurf_t0000000.vtk UpperSurf_n0.vtk
yes | cp -rf UpperSurf_t0000000.vtk $cwd/Results/UpperSurf_n0s.vtk
##Get Strain Profile and Create new geometry
yes | cp -rf UpperSurf_t0000000.vtk UpperSurf_thick.vtk
python Process_Thick.py
yes | cp -rf Strain_Profile.txt $cwd
yes | cp -rf Strain_Profile.txt $cwd/Results/Strain_Profile_n0s.txt
cd $cwd




value=1
count=100
for i in $(seq $count); do

	yes | cp -rf geometry_base.exo geometry.exo
	python Modify_Thickness_Strain_Prev.py
	yes | cp -rf Strain_Profile.txt Strain_Profile_Prev.txt
	yes | cp -rf Thickness_Profile.txt $cwd/Results/Thickness_Profile_n"$((i-1))".txt


	source Run_Pylith_dummy.sh
	source shell/fix_path.sh
	cd output 
	yes | cp -rf All_Nodes_t0000000.vtk $cwd
	yes | cp -rf All_Nodes_t0000000.vtk $cwd/Results/All_Nodes_n"$((i-1))".vtk
	yes | cp -rf LeftSurf_t0000000.vtk $cwd
	yes | cp -rf RightSurf_t0000000.vtk $cwd
	cd $cwd
	python Make_ForceFile.py
	yes | cp -rf forces.spatialdb spatialdb/forces.spatialdb
	source Run_Pylith.sh
	source shell/fix_path.sh
	cd output
	yes | cp -rf UpperSurf_t0000000.vtk $cwd
	yes | cp -rf UpperSurf_t0000000.vtk UpperSurf_n0.vtk
	yes | cp -rf UpperSurf_t0000000.vtk $cwd/Results/UpperSurf_n"$((i-1))".vtk
	##Get Strain Profile and Create new geometry
	yes | cp -rf UpperSurf_t0000000.vtk UpperSurf_thick.vtk
	python Process_Thick.py
	yes | cp -rf Strain_Profile.txt $cwd
	yes | cp -rf Strain_Profile.txt $cwd/Results/Strain_Profile_n"$((i-1))".txt
	cd $cwd
	yes | cp -rf geometry_base.exo geometry.exo
	python Modify_Thickness_Strain_Prev.py
	yes | cp -rf Strain_Profile.txt Strain_Profile_Prev.txt


done



##

cd $cwd/Results
python Plot_Results.py








