# Workflow Example 4: Irregular domain with a river mask
# This example walks through a case where
# you have an irregular domain and a pre-defined river network to 
# be used for the DEM processing 
# This requires two inputs: (1) a DEM  (2) a mask of the domain and (3) a river mask

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
# Inputs
# The DEM and mask should be formated as a matrices with the same 
# dimensions as the domain (i.e. ncol=nx, nrow=ny)
# The river and domain masks should consist of 0's and 1's with 1's for any grid cell on the river network or inside the domain respectively
dem=matrix(scan("dem_test.txt"), ncol=215, byrow=T)
rivermask=matrix(scan("river_mask_test.txt"), ncol=215, byrow=T) #mask of river cells
domainmask=matrix(scan("mask_test.txt"), ncol=215, byrow=T) #Mask of domain extent

ny=nrow(dem)
nx=ncol(dem)

#transforming the inputs so it is indexed as [x,y] for functions
demT=t(dem[ny:1,])
rivermaskT=t(rivermask[ny:1,])
domainmaskT=t(domainmask[ny:1,])
rivermaskT=rivermaskT*domainmaskT #clipping the river mask to the domain mask

#check that you aren't upside down and backwards somehow...
#if you've formatted your input correctly this should look 
#like your domain without any additional transforming
image.plot(demT)
image.plot(rivermaskT)
image.plot(domainmaskT)

##########################################
# Process the DEM: 
#1. Initialize the queue with river cells that fall on the border
#2. Traverse the stream network filling sinks and stair stepping around D8 neigbhors
#3. Look for orphan branches and continue processing until they are all connected
#4. Use the processed river cells as the intialize a new queue 
#5. process hillslopes from there
ep=0.01 #epsilon to add to flat cell

#1.initialize the queue with river cells that fall on the border
init=InitQueue(demT,  initmask=rivermaskT, domainmask=domainmaskT) 

#2.take a first pass at traversing the streams
trav1=StreamTraverse(demT, mask=rivermaskT, init$queue, init$marked, basins=init$basins, printstep=F, epsilon=ep)
print(paste("First Pass:", round(sum(trav1$marked * rivermaskT)/sum(rivermaskT)*100,1), " % cells processed"))

image(trav1$basins, zlim=c(0.5, max(trav1$basins)))
image.plot(trav1$marked) #The portion of the river mask traversed so far


#3. look for 'orphaned' branches and continue traversing until they are all connected
# orphaned branches are portions of the river network that are connected diagonally (i.e. without any d4 neighbors)
norphan=1
lap=1
while(norphan>0){
	#look for orphan branches
	orphan=FindOrphan(trav1$dem, rivermaskT, trav1$marked)
	norphan=orphan$norphan
	print(paste("lap", lap, norphan, "orphans found"))
	
	#go around again if orphans are found
	if(norphan>0){
		trav2 = StreamTraverse(trav1$dem, mask=rivermaskT, queue=orphan$queue, marked=trav1$marked, step=trav1$step, direction=trav1$direction, basins=trav1$basins, printstep=F, epsilon=ep)
		trav1=trav2
		lap=lap+1
	} else {
		print("Done!  No orphan branches found")
	}
} #end while
print(paste("Final pass:", round(sum(trav1$marked * rivermaskT)/sum(rivermaskT)*100,1), " % cells processed"))

image(trav1$mask)
image(trav1$basins, zlim=c(0.5, max(trav1$basins)))
image(trav1$marked)


#4.initialize the queue with every cell on the processed river boundary
#to do this use the marked rivers from the last step plus the edge cells
#as the boundary and the mask
borderT=GetBorder(domainmaskT)
RivBorder=borderT+trav1$marked 
RivBorder[RivBorder>1]=1
image(RivBorder)
init=InitQueue(trav1$dem,  border=RivBorder)


#5.process all the cells off the river usins the river as the boundary
travHS=D4TraverseB(trav1$dem, init$queue, init$marked, direction=trav1$direction, basins=trav1$basins, step=trav1$step, epsilon=ep, mask=domainmaskT) #irregular boundary

image(travHS$marked) #this shoudl be the entire mask now
image(travHS$step)
image(travHS$basins,zlim=c(0.5, max(travHS$basins)), col=grey.colors(10))
maskcol=colorRampPalette(c('white', 'red'))
image(trav1$marked,zlim=c(0.5,1), col=maskcol(2), add=T)

##########################################
# Calculate the slopes
# Note this step also fixes the directions of the borders because
# directions are not provided when the queue is initialized

# Option 1: just calcualte the slopes for the entire domain with no distinction between river and hillslope cells
#In this example secondary slope scaling is turned on and the secondary
#Slopes are set to a maximum of 0.1*primary flow direction
#To calculate only slopes in the primary flow direction set the secondaryTH to 0
slopesUW=SlopeCalcUP(dem=travHS$dem, direction=travHS$direction, dx=1000, dy=1000, scalesecond=T, secondaryTH=0.1, mask=domainmaskT) 

#Option 2: If you would like to scale the secondary slopes differently for hillslopes than for rivers 
#Make a mask of cells with drainage area greater than some threshold to define as the rivers (you could also use the original river mask here)
area=drainageArea(travHS$direction, mask=domainmaskT, printflag=T) #rectangular boundary
image.plot(area)
rth=50 #Drainage area threshold in number of grid cells
rivers=area
rivers[area<rth]=0
rivers[area>=rth]=1
image.plot(rivers) #Plot to check that the threshold is good

slopesUW=SlopeCalcUP(dem=travHS$dem, mask=domainmaskT, direction=travHS$direction, dx=1000, dy=1000, scalesecond=T, secondaryTH=0.1, rivermask=rivers) 


#Look at the slopes and directions
image(slopesUW$slopex)
image(slopesUW$slopey)
image.plot(slopesUW$direction)


##########################################
# Calculate the drainage area
area=drainageArea(slopesUW$direction, mask=domainmaskT, printflag=T) 
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



