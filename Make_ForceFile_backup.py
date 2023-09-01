

##File to extract node information and apply boundary forces....

import numpy as np
import pandas as pd
from utility import Sph2Cart_2
from utility import Cart2Sph_2
from utility import Sph2Cart
from utility import Cart2Sph
from utility import Get_Cellvols_2D
from utility import Get_Cellvols_3D
from utility import Parse_vtk_2D
from utility import Parse_vtk_3D
from utility import Get_NodeVols
import pickle



finame="All_Nodes_t0000000"
NodesData,CellsData=Parse_vtk_2D.main(finame)[0],Parse_vtk_2D.main(finame)[1]

print("Forces to Sides...")

Cellvol=Get_Cellvols_2D.main(CellsData,NodesData)[0]
TotalVol=np.sum(Cellvol)
NodeVol=Get_NodeVols.main(NodesData,CellsData,Cellvol)


xpoints=NodesData[:,0]
ypoints=NodesData[:,1]
zpoints=NodesData[:,2]

veclen=len(xpoints)
        
Traction= -10000*xpoints
TotNodeVol=np.sum(NodeVol)
NodeVolReal=((NodeVol)/TotNodeVol)*TotalVol

Fin_Traction=np.transpose(NodeVolReal) * Traction  ##Stress times area to get Force



###Need to zero out values for nodes not on the edges

Edge_Length=np.loadtxt("Edge_Length.txt")

for i in range(veclen):
	if xpoints[i]>-Edge_Length/2+1 and xpoints[i]<Edge_Length/2-1:
		Fin_Traction[:,i] = 0
		#print(Fin_Traction[:,i])
	#print(xpoints[i])
	#print(-Edge_Length/2+1)
		





Forcefun=np.zeros([veclen,4])

for x in range(veclen):
	Forcefun[x,0]=NodesData[x,0];  ##keep everything in m
	Forcefun[x,1]=NodesData[x,1];
	Forcefun[x,2]=Fin_Traction[:,x] 
	Forcefun[x,3]=0;


##Outputting

SpatialFile=Forcefun;


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




