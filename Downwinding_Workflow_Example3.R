# Downwinding Workflow Example 3: 
# - A rectangular domain and a pre-defined river network to be used for the DEM processing
# - Calculating slopes with a downwinding approach to be consistent with ParFlow's OverlandFlow boundary conditon. 
# - This requires two inputs: (1) a DEM and (2) a river mask

#  This example usese the test domain from Condon and Maxwell (2019) (https://doi.org/10.1016/j.cageo.2019.01.020)
#  the datasets for this domain are provided with the PriorityFlow R package
#  to use your own datasets you should have a DEM and mask files formatted as a matrices
#  with [i,j] corresponding to the x and y axes of the domain 
#  (i.e. DEM[1,1] is the lower left corner of the domain and DEM[nx,ny] is the upper right)
#  See help(DEM) and help(river.mask) for more information

##########################################
#Source Libraries and functions
#Refer to the ReadMe on the Github repo for how to install the 
#PriorityFlow library if you haven't installed it yet: https://github.com/lecondon/PriorityFlow
rm(list=ls())
library("PriorityFlow")
help("PriorityFlow")
library('fields')

##########################################
#Settings
#Settings for the slope processing to change

#DEM processing
ep=0.01 #The epsilon value applied to flat cells

#Slope scaling
maxslope=0.5	#maximum slope (slopes with absolute value greater than this will be set to minslope), set to -1 if you don't want to set a threshold
minslope=1e-5	#minimum slope (slopes with absolute value less than this will be set to minslope), set to -1 if you don't want to set a threshold
scale=0.1 #The maximum ratio of secondary to primary flow directions (set to -1 if you don't want the secondary slopes to be scaled, set to 0 if you want only primary flow directios)

#River and subbasin size for slope calculations
sub_th=100 #area threshold (number of grid cells) to use for subbasin delineation
riv_th=sub_th #optional additional area threshold (number of grid cells) to use for the river mask for slope processing. See notes below if you want to change this to be different from the subbasin threshold
riv_method=3 #method for processing river cellls (0=treat river cells the same as the rest of the domain, 1=set secondary slopes along the river to zero, 2=apply subbasin average slopes to river cells, 3=apply subbasin average river slopes of the river cells)
mrg_th=10	#Threshold number of grid cells for merging small subbasins in the subbasin analysis


#Grid dimensions for slopes
dx=1000 #grid cell size for slope calculations
dy=1000 #grid cell size for slope calcualtions

#Runname - for output writing
runname='Test'

##########################################
# Get the dimensions of the domain and look at the inputs
nx=nrow(DEM)
ny=ncol(DEM)

par(mfrow=c(1,2))
image.plot(DEM)
image.plot(river.mask)

##########################################
# Process the DEM:
#1. Initialize the queue withriver cells that fall on the border
#2. Traverse the stream network filling sinks and stair stepping around D8 neigbhors
#3. Look for orphan branches and continue processing until they are all connected
#4. Use the processed river cells as the intialize a new queue
#5. process hillslopes from there

#1.initialize the queue with river cells that fall on the border
init=InitQueue(DEM,  initmask=river.mask) #rectangular boundary

#2.take a first pass at traversing the streams
trav1=StreamTraverse(DEM, mask=river.mask, init$queue, init$marked, basins=init$basins, printstep=F, epsilon=ep)
print(paste("First Pass:", round(sum(trav1$marked * river.mask)/sum(river.mask)*100,1), " % cells processed"))

image(trav1$basins, zlim=c(0.5, max(trav1$basins)))
image.plot(trav1$marked) #The portion of the river mask traversed so far


#3. look for 'orphaned' branches and continue traversing until they are all connected
# orphaned branches are portions of the river network that are connected diagonally (i.e. without any d4 neighbors)
norphan=1
lap=1
while(norphan>0){
	#look for orphan branches
	orphan=FindOrphan(trav1$dem, river.mask, trav1$marked)
	norphan=orphan$norphan
	print(paste("lap", lap, norphan, "orphans found"))

	#go around again if orphans are found
	if(norphan>0){
		trav2 = StreamTraverse(trav1$dem, mask=river.mask, queue=orphan$queue, marked=trav1$marked, step=trav1$step, direction=trav1$direction, basins=trav1$basins, printstep=F, epsilon=ep)
		trav1=trav2
		lap=lap+1
	} else {
		print("Done!  No orphan branches found")
	}
} #end while
print(paste("Final pass:", round(sum(trav1$marked * river.mask)/sum(river.mask)*100,1), " % cells processed"))

image(trav1$mask)
image(trav1$basins, zlim=c(0.5, max(trav1$basins)))
image(trav1$marked)


#4.initialize the queue with every cell on the processed river boundary
#to do this use the marked rivers from the last step plus the edge cells
#as the boundary and the mask
inittemp=InitQueue(DEM) # Rectangular domain just using this to get a marked matrix of the edge cells
RivBorder=inittemp$marked+trav1$marked #initializing with rectangular boundary
RivBorder[RivBorder>1]=1
image(RivBorder)
init=InitQueue(trav1$dem,  border=RivBorder)

#5.process all the cells off the river usins the river as the boundary
travHS=D4TraverseB(trav1$dem, init$queue, init$marked, direction=trav1$direction, basins=trav1$basins, step=trav1$step, epsilon=ep)

image(travHS$marked) #this should be the entire domain now
image(travHS$step)
image(travHS$basins,zlim=c(0.5, max(travHS$basins)), col=grey.colors(10))
maskcol=colorRampPalette(c('white', 'red'))
image(trav1$marked,zlim=c(0.5,1), col=maskcol(2), add=T)

##########################################
# Calculate the slopes
# Note this step also fixes the directions of the borders because
# directions are not provided when the queue is initialized

### Option 1: just calcualte the slopes for the entire domain with no distinction between river and hillslope cells
#In this example secondary slope scaling is turned on and the secondary
#Slopes in the secondary direction are set to a maximum of 0.1*primary flow direction
#To calculate only slopes in the primary flow direction set the secondaryTH to 0
#Additionally primary slopes are limited by min slope and max slope thresholds
slopesUW=SlopeCalcUP(dem=travHS$dem, direction=travHS$direction, dx=dx, dy=dy,  secondaryTH=scale, maxslope=maxslope, minslope=minslope)

### Option 2: If you would like to handle river cells differently from the rest of the domain-- NOTE the 'river cells' here will be determined based on drainage area thresholds after not the orinigal rive mask provided.
# Calculate the drainage area
area=drainageArea(travHS$direction, printflag=F)

# Define subbasins for calcualting river reach slopes
# the riv_th here is the drainage area threshold for splitting the river network branches
# when you do this you can still end up with subbasins with drainage areas less than the riv_th
# when multiple branches come together in a small area.
# To fix this you can set a merge threshold (merge_th) so that subbains with areas < merge_th autmoatically get merged with their downstream neighbor
subbasin=CalcSubbasins(travHS$direction, area, riv_th=sub_th, merge_th=mrg_th)
#plot the resulting subbasins and rivers
temp=subbasin$RiverMask
temp[temp==0]=NA
maskcol=colorRampPalette(c('black', 'black'))
#maskcol=colorRampPalette(c('white', 'white'))
image.plot(subbasin$subbasins)
image.plot((temp*2), add=T, col=maskcol(2), legend=F)

#Calculate the slopes
# The "river_method' flag here determines how the river cells will be handeled (e.g. using subbasin averages along reaches). Refer to the top of this script or the function for details.
slopesUW=SlopeCalcUP(dem=travHS$dem, direction=travHS$direction, dx=dx, dy=dy, secondaryTH=scale, maxslope=maxslope, minslope=minslope, river_method=riv_method, rivermask=subbasin$RiverMask, subbasin=subbasin$subbasins)

#Alternate more advanced approach: Define a river mask separate from the subbasin river mask and use this for the slope calculations. If you do this the average slopes will still be calculated
#by subbasin using the sub_th, but you can apply those average sloeps to more river cells by setting a lower threshold here. This is the 'riv_th' set at the top
#if you set riv_th=sub_th at the top this will have the same effect as just running the slope calc with the subbasin$RiverMask
rivers=area
rivers[area<riv_th]=0
rivers[area>=riv_th]=1

#plot the subbasins with the new river mask to check that the threshold is good
temp=rivers
temp[temp==0]=NA
maskcol=colorRampPalette(c('black', 'black'))
#maskcol=colorRampPalette(c('white', 'white'))
image.plot(subbasin$subbasins)
image.plot((temp*2), add=T, col=maskcol(2), legend=F)

slopesUW=SlopeCalcUP(dem=travHS$dem, direction=travHS$direction, dx=dx, dy=dy, secondaryTH=scale, maxslope=maxslope, minslope=minslope, river_method=riv_method, rivermask=rivers, subbasin=subbasin$subbasins)

#Look at the slopes and directions
image(slopesUW$slopex)
image(slopesUW$slopey)
image.plot(slopesUW$direction)

##########################################
# Calculate the drainage area
area=drainageArea(slopesUW$direction, printflag=F)
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

fout=paste(runname, ".slopex.sa", sep="")
write.table( t(c(nx,ny,1)), fout, append=F, row.names=F, col.names=F)
write.table(slopeUWx, fout, append=T, row.names=F, col.names=F)
fout=paste(runname,".slopey.sa", sep="")
write.table( t(c(nx,ny,1)), fout, append=F, row.names=F, col.names=F)
write.table(slopeUWy, fout, append=T, row.names=F, col.names=F)


##########################################
#Example writing out other variables as matrices
write.table( t(slopesUW$direction[,ny:1]) ,paste(runname, ".direction.out.txt", sep=""), row.names=F, col.names=F)
write.table( t(travHS$dem[,ny:1]) ,paste(runname, ".dem.out.txt", sep=""), row.names=F, col.names=F)
write.table( t(area[,ny:1]) , paste(runname, ".area.out.txt", sep=""), row.names=F, col.names=F)
write.table( t(subbasin$subbasins[,ny:1]) , paste(runname, ".subbasins.out.txt", sep=""), row.names=F, col.names=F)
write.table( t(subbasin$segments[,ny:1]) , paste(runname, ".subbasin_streams.out.txt", sep=""), row.names=F, col.names=F)

## write out the subbasin summary information
write.table(subbasin$summary, paste(runname, ".Subbasin_Summary.txt", sep=""), row.names=F)
