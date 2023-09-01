

# *********************************************************************
# FUNCTION TO EVALUATE CELL AREAS (2D) FROM MESH INFORMATION
#
# *********************************************************************

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
import pandas as pd




def main(CellsData,NodesData):


	

	Cellsize_arr=len(CellsData[:,1]);
	Cellloc=np.zeros([Cellsize_arr,3]); ##Three Sptial Dimensions
	Cellvol=np.zeros([Cellsize_arr,1]);
	Nodesize_arr=len(NodesData[:,2]);
	Cellvol2=np.zeros([Cellsize_arr,1]);

	for x in range(Cellsize_arr):
		Node1num= CellsData[x,0];
		Node2num= CellsData[x,1];
		Node1xval=NodesData[Node1num,0];
		Node1yval=NodesData[Node1num,1];
		Node2xval=NodesData[Node2num,0];
		Node2yval=NodesData[Node2num,1];
		xloc=(Node1xval+Node2xval)/2 ; 
		yloc=(Node1yval+Node2yval)/2 ;
		##Calculate Location
		Cellloc[x,0]=xloc;
		Cellloc[x,1]=yloc;
		##Calculate Area in m of each triangle
		Cellvol[x]=((Node1xval-Node2xval)**2+(Node1yval-Node2yval)**2)**(1/2)


	return Cellvol,Cellloc


