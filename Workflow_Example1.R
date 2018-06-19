# Workflow Example 1: Rectangular domain no river mask
# This example walks through the simplest case where
# you have a rectangular domain and no river network specified a-priori
# in this case the only input required is a DEM

##########################################
#Source Libraries and functions
rm(list=ls())
library('fields')
source("./functions/D4_Traverse.R")
source("./functions/Init_Queue.R")
source("./functions/Stream_Traverse.R")
source("./functions/Find_Orphan.R")
source("./functions/drainage_area.R")
source("./functions/Slope_Calc_Upwind.R")
source("./functions/Get_Border.R")

##########################################
# Read in Inputs
# The DEM should be formated as a matrix with the same 
# dimensions as the domain
dem=matrix(scan("dem_test.txt"), ncol=215, byrow=T)
ny=nrow(dem)
nx=ncol(dem)
demT=t(dem[ny:1,] #transforming the dem so it is indexed as [x,y] for functions

#check that you aren't upside down and backwards somehow...
#if you've formatted your input correctly this should look like your DEM
# without any additional transforming on the matrix
image.plot(demT)

##########################################
# Process the DEM
#1. initialize the queue with all the rectangular border cells
#2. Process the DEM so that all cells drain to the boundaries

#1. Initialize queue
init=InitQueue(demT) #using the rectangular boundary
#2. Process DEM
travHS=D4TraverseB(demT, init$queue, init$marked, basins=init$basins, epsilon=0.01) 

#Look at the outputs
ls(travHS) #to see a list of everything that comes out of the processing
image.plot(travHS$dem) # the processed DEM
image(travHS$marked) # Mask of cells the processing algorithm covered (should be everything for this approach)
image(travHS$step) # The step at which each cell was processed in the algorithm
image(travHS$basins) # The resulting drainage basins

##########################################
# Calculate the slopes
# Note this step also fixes the directions of the borders because
# directions are not provided when the queue is initialized

# Option 1: just calcualte the slopes for the entire domain with no distinction between river and hillslope cells
#In this example secondary slope scaling is turned on and the secondary
#Slopes are set to a maximum of 0.1*primary flow direction
#To calculate only slopes in the primary flow direction set the secondaryTH to 0
slopesUW=SlopeCalcUP(dem=travHS$dem, direction=travHS$direction, dx=1000, dy=1000, scalesecond=T, secondaryTH=0.1) 


#Option 2: If you would like to scale the secondary slopes differently for 
#Hillslopes than for rivers 
#Make a mask of cells with drainage area greater than some threshold to define as the rivers (you could also input your own river mask here)
area=drainageArea(travHS$direction, printflag=T) #rectangular boundary
image.plot(area)
rth=50 #Drainage area threshold in number of grid cells
rivers=area
rivers[area<rth]=0
rivers[area>=rth]=1
image.plot(rivers) #Plot to check that the threshold is good

slopesUW=SlopeCalcUP(dem=travHS$dem, direction=travHS$direction, dx=1000, dy=1000, borderdir=1, scalesecond=T, secondaryTH=0.1, rivermask=rivers) 


#Look at the slopes and directions
image(slopesUW$slopex)
image(slopesUW$slopey)
image.plot(slopesUW$direction)

##########################################
# Calculate the drainage area
area=drainageArea(slopesUW$direction, printflag=T) #rectangular boundary
image.plot(area)

##########################################
#Write the slopes out in PF format
slopeUWx=slopeUWy=rep(0, nx*ny)
jj=1
for(j in 1:ny){
	for(i in 1:nx){
		slopeUWx[jj]=slopesUW$slopex[i,j]
		slopeUWy[jj]=slopesUW$slopey[i,j]
		jj=jj+1
	}
}

fout="slopex.sa"
write.table( t(c(nx,ny,1)), fout, append=F, row.names=F, col.names=F)
write.table(slopeUWx, fout, append=T, row.names=F, col.names=F)
fout="slopey.sa"
write.table( t(c(nx,ny,1)), fout, append=F, row.names=F, col.names=F)
write.table(slopeUWy, fout, append=T, row.names=F, col.names=F)

##########################################
#Example writing out other variables as matrices 
write.table( t(travHS$direction[,ny:1]) ,"direction.out.test.txt", row.names=F, col.names=F)
write.table( t(travHS$dem[,ny:1]) ,"dem.out.test.txt", row.names=F, col.names=F)
write.table( t(area[,ny:1]) ,"area.out.test.txt", row.names=F, col.names=F)