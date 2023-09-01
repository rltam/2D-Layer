

##File to extract node information and apply boundary forces....

import numpy as np
import pandas as pd
from utility import Sph2Cart_2
from utility import Cart2Sph_2
from utility import Sph2Cart
from utility import Cart2Sph
from utility import Get_Cellvols_1D
from utility import Get_Cellvols_2D
from utility import Get_Cellvols_3D
from utility import Parse_vtk_1D
from utility import Parse_vtk_2D
from utility import Parse_vtk_3D
from utility import Get_NodeVols
from utility import Get_NodeVols_1D
import pickle



finame="LeftSurf_t0000000"
NodesData,CellsData=Parse_vtk_1D.main(finame)[0],Parse_vtk_1D.main(finame)[1]

NodesData1=NodesData
print("Forces to Sides...")

Cellvol=Get_Cellvols_1D.main(CellsData,NodesData)[0]
TotalVol=np.sum(Cellvol)
NodeVol=Get_NodeVols_1D.main(NodesData,CellsData,Cellvol)


xpoints=NodesData[:,0]
ypoints=NodesData[:,1]
zpoints=NodesData[:,2]

veclen=len(xpoints)
      
Traction= 5000*xpoints
TotNodeVol=np.sum(NodeVol)
NodeVolReal=((NodeVol)/TotNodeVol)*TotalVol
Fin_Traction1=np.transpose(NodeVolReal) * Traction  ##Stress times area to get Force


###Need to zero out values for nodes not on the edges

Edge_Length=np.loadtxt("Edge_Length.txt")

#for i in range(veclen):
#	if xpoints[i]>-Edge_Length/2+1 and xpoints[i]<Edge_Length/2-1:
#		Fin_Traction[:,i] = 0
		#print(Fin_Traction[:,i])
	#print(xpoints[i])
	#print(-Edge_Length/2+1)
		

finame="RightSurf_t0000000"
NodesData,CellsData=Parse_vtk_1D.main(finame)[0],Parse_vtk_1D.main(finame)[1]

NodesData2=NodesData
print("Forces to Sides...")

Cellvol=Get_Cellvols_1D.main(CellsData,NodesData)[0]
TotalVol=np.sum(Cellvol)
NodeVol=Get_NodeVols_1D.main(NodesData,CellsData,Cellvol)


xpoints=NodesData[:,0]
ypoints=NodesData[:,1]
zpoints=NodesData[:,2]

veclen=len(xpoints)
      
Traction= 10000*xpoints*1e4
TotNodeVol=np.sum(NodeVol)
NodeVolReal=((NodeVol)/TotNodeVol)#*TotalVol
Fin_Traction2=np.transpose(NodeVolReal) * Traction  ##Stress times area to get Force



Fin_Traction=np.concatenate((Fin_Traction1,Fin_Traction2),axis=1)
NodesData=np.concatenate((NodesData1,NodesData2),axis=0)


veclen=len(NodesData[:,0])

Forcefun=np.zeros((veclen,4))

for x in range(veclen):
	Forcefun[x,0]=NodesData[x,0]  ##keep everything in m
	Forcefun[x,1]=NodesData[x,1]
	Forcefun[x,2]=Fin_Traction[:,x] 
	Forcefun[x,3]=0



###We need to assign zero to all other values....

Zero_trac=np.zeros((len(Forcefun[:,0]),4))

print(Zero_trac.shape)


for i in range(len(Zero_trac)):
	if NodesData[i,0]<0:
		Zero_trac[i,0]=NodesData[i,0]+1
		Zero_trac[i,1]=NodesData[i,1]
		Zero_trac[i,2]=0
		Zero_trac[i,3]=0
	elif NodesData[i,0]>0:
		Zero_trac[i,0]=NodesData[i,0]-1
		Zero_trac[i,1]=NodesData[i,1]
		Zero_trac[i,2]=0
		Zero_trac[i,3]=0

veclen=len(Zero_trac[:,0])


Forcefun2=np.zeros([veclen,4])

for x in range(veclen):
	Forcefun2[x,0]=Zero_trac[x,0]  ##keep everything in m
	Forcefun2[x,1]=Zero_trac[x,1]
	Forcefun2[x,2]=Zero_trac[x,2] 
	Forcefun2[x,3]=Zero_trac[x,3]



##Outputting


ForceFun=np.concatenate((Forcefun,Forcefun2),axis=0)

veclen=len(ForceFun[:,0])



SpatialFile=ForceFun;


Spatialdb=open("forces.spatialdb","w")
Spatialdb.write("#SPATIAL.ascii 1\n")
Spatialdb.write("SimpleDB {\n")
Spatialdb.write("  num-values = 2\n")
Spatialdb.write("  value-names = force-x force-y\n")
Spatialdb.write("  value-units = m m\n")
Spatialdb.write("  num-locs = ")
Spatialdb.write(str(veclen))
Spatialdb.write("\n")
Spatialdb.write("  data-dim = 2\n")
Spatialdb.write("  space-dim = 2\n")
Spatialdb.write("  cs-data = cartesian {\n")
Spatialdb.write("    to-meters = 1.0e+0\n")
Spatialdb.write("    space-dim = 2\n")
Spatialdb.write("  }\n")
Spatialdb.write("}\n")


np.savetxt(Spatialdb,SpatialFile,fmt="%d")





