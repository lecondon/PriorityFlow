# Downwinding Workflow Example 2: 
# - Irregular domain boundary and no river mask
# - Calculating slopes with a downwinding approach to be consistent with ParFlow's OverlandFlow boundary conditon. 
# - This requires two inputs: (1) a DEM and (2) a watershed mask

# This example usese the test domain from Condon and Maxwell (2019) (https://doi.org/10.1016/j.cageo.2019.01.020)
#  the datasets for this domain are provided with the PriorityFlow R package
#  to use your own datasets you should have a DEM and mask files formatted as a matrices
#  with [i,j] corresponding to the x and y axes of the domain 
#  (i.e. DEM[1,1] is the lower left corner of the domain and DEM[nx,ny] is the upper right)
#  See help(DEM) and help(watershed.mask) for more information

##########################################
#Source Libraries and functions
#Note refer to the ReadMe on the Github repo for how to install the 
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
riv_th=sub_th #optional additional area threshold (number of grid cells)  to use for the river mask for slope processing. See notes below if you want to change this to be different from the subbasin threshold
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
image.plot(watershed.mask)


##########################################
# Process the DEM:
#1. initialize the queue with all the border cells for the irrgular boundary
#2. Process the DEM so that all cells drain to the borders

#1. Initialize queue
init=InitQueue(DEM, domainmask=watershed.mask) #using the rectangular boundary
#2. Process DEM
travHS=D4TraverseB(DEM, init$queue, init$marked, basins=init$basins, epsilon=ep, mask=watershed.mask)

#Look at the outputs
ls(travHS) #to see a list of everything that comes out of the processing
image.plot(travHS$dem) # the processed DEM
image(travHS$marked) # Mask of cells the processing algorithm covered (should be everything inside the mask)
image(travHS$step) # The step at which each cell was processed in the algorithm
image(travHS$basins) # The resulting drainage basins


##########################################
# Calculate the slopes
# Note this step also fixes the directions of the borders because
# directions are not provided when the queue is initialized

### Option 1: just calcualte the slopes for the entire domain with no distinction between river and hillslope cells
#In this example secondary slope scaling is turned on and the secondary
#Slopes in the secondary direction are set to a maximum of 0.1*primary flow direction
#To calculate only slopes in the primary flow direction set the secondaryTH to 0
#Additionally primary slopes are limited by min slope and max slope thresholds
slopesUW=SlopeCalcUP(dem=travHS$dem, mask=watershed.mask, direction=travHS$direction, dx=dx, dy=dy,  secondaryTH=scale, maxslope=maxslope, minslope=minslope)


### Option 2: If you would like to handle river cells differently from the rest of the domain
# Calculate the drainage area
area=drainageArea(travHS$direction, printflag=F)

# Define subbasins for calcualting river reach slopes
# the riv_th here is the drainage area threshold for splitting the river network branches
# when you do this you can still end up with subbasins with drainage areas less than the riv_th
# when multiple branches come together in a small area.
# To fix this you can set a merge threshold (merge_th) so that subbains with areas < merge_th autmoatically get merged with their downstream neighbor
subbasin=CalcSubbasins(travHS$direction, mask=watershed.mask, area, riv_th=sub_th, merge_th=mrg_th)
#plot the resulting subbasins and rivers
temp=subbasin$RiverMask
temp[temp==0]=NA
maskcol=colorRampPalette(c('black', 'black'))
#maskcol=colorRampPalette(c('white', 'white'))
image.plot(subbasin$subbasins)
image.plot((temp*2), add=T, col=maskcol(2), legend=F)

#Calculate the slopes
# The "river_method' flag here determines how the river cells will be handeled (e.g. using subbasin averages along reaches). Refer to the top of this script or the function for details.
slopesUW=SlopeCalcUP(dem=travHS$dem, mask=watershed.mask, direction=travHS$direction, dx=dx, dy=dy, secondaryTH=scale, maxslope=maxslope, minslope=minslope, river_method=riv_method, rivermask=subbasin$RiverMask, subbasin=subbasin$subbasins)

### Option 2b: Alternate more advanced approach: Define a river mask separate from the subbasin river mask and use this for the slope calculations. If you do this the average slopes will still be calculated
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

slopesUW=SlopeCalcUP(dem=travHS$dem,  mask=watershed.mask, direction=travHS$direction, dx=dx, dy=dy, secondaryTH=0.1, maxslope=maxslope, minslope=minslope, river_method=riv_method, rivermask=rivers, subbasin=subbasin$subbasins)


#Look at the slopes and directions
image(slopesUW$slopex)
image(slopesUW$slopey)
image.plot(slopesUW$direction)


##########################################
# Calculate the drainage area - if you went with option 1 for slopes and you didn't do this already
area=drainageArea(travHS$direction, mask=watershed.mask, printflag=F)
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
write.table( t(travHS$direction[,ny:1]) ,paste(runname, ".direction.out.txt", sep=""), row.names=F, col.names=F)
write.table( t(travHS$dem[,ny:1]) ,paste(runname, ".dem.out.txt", sep=""), row.names=F, col.names=F)
write.table( t(area[,ny:1]) , paste(runname, ".area.out.txt", sep=""), row.names=F, col.names=F)
write.table( t(subbasin$subbasins[,ny:1]) , paste(runname, ".subbasins.out.txt", sep=""), row.names=F, col.names=F)
write.table( t(subbasin$segments[,ny:1]) , paste(runname, ".subbasin_streams.out.txt", sep=""), row.names=F, col.names=F)

## write out the subbasin summary information
write.table(subbasin$summary, paste(runname, ".Subbasin_Summary.txt", sep=""), row.names=F)
