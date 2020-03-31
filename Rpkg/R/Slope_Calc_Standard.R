#'Calculate the slopes from a DEM
#'
#'This function will calculate slopes using standard or upwinding options and apply a range of smoothing options

#' @inheritParams D4TraverseB
#' @param  dx,dy Lateral grid cell resolution
#' @param minslope Minimum absolute slope value to apply to flat cells if needed. Defaults to 1e-5
#' @param maxslope  Maximum absolute value of slopes. If this is set to -1 the slopes will not be limited. Default value is -1
#' @param secondaryTH  secondary threshold - maximum ratio of |secondary|/|primary| to be enforced. 
#' NOTE - this scaling occurs after any max threholds are applied. Currently this is only working for two options:
#' (1) If this is set to -1 no scaling will be applied, (2) If  this is set to zero all  seconeary slopes will be zero
#' @param printflag Print function progress

#' @section River Methods:
#' 0: default value, no special treatment for river cells
#'
#' 1: Scale secondary slopes along the river (Note this requries a river mask and you must set a river_secondaryTH if you want this to be something other than 0)
#'
#' 2: Apply watershed mean slope to each river reach (requires river mask and subbasins)
#'
#' 3: Apply the stream mean slope to each reach (requires river mask and subbasins)
#'
#' NOTE: the river mask can be different from the rivers that were used to create the subbasins if desired (i.e. if you want to use a threshold of 100 to create subbasins but then apply to river cells with a threshold of 50)
#' @export
SlopeCalStan=function(dem, direction, dx, dy,  mask,  d4=c(1,2,3,4),  minslope=0, maxslope=-1, secondaryTH=-1, printflag=F){
####################################################################
# PriorityFlow - Topographic Processing Toolkit for Hydrologic Models
# Copyright (C) 2018  Laura Condon (lecondon@email.arizona.edu)
# Contributors - Reed Maxwell (rmaxwell@mines.edu)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation version 3 of the License

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>
####################################################################


###############################################################
ny=ncol(dem)
nx=nrow(dem)

#If no mask is provided default to the rectugular domain
if(missing(mask)){
	print("No mask provided, initializing mask for complete domain")
	mask=matrix(1, nrow=nx, ncol=ny)
	borders=matrix(0, nrow=nx, ncol=ny)
}

#Identify the border cells
borders=matrix(1, nrow=nx, ncol=ny)
borders[2:(nx-1), 2:(ny-1)]= mask[1:(nx-2), 2:(ny-1)] + 																
	              mask[3:nx, 2:(ny-1)] +
								mask[2:(nx-1), 1:(ny-2)] +
								mask[2:(nx-1), 3:ny]
borders=borders*mask
borders[which(borders<4 & borders!=0)]=1 
borders[borders==4]=0


bordi=borders #making a border indicator file with all 1's for any border cell
bordi[bordi>1]=1
bordlist=which(borders>0)

#assign NA values to the DEM outside the mask
inmask=which(mask==1)
demMask=dem
demMask[-inmask]=NA


#First pass calculate the x and y slopes as (i+1)-(i)
# slopex = dem[i+1,j]- dem[i,j]
# slopey = dem[i,j+1] - dem[i,j]
slopex1=matrix(NA, ncol=ny, nrow=nx)
slopey1=matrix(NA, ncol=ny, nrow=nx)
slopex1[1:(nx-1),]=(demMask[2:nx,]-demMask[1:(nx-1),]) /dx
slopey1[,1:(ny-1)]=(demMask[,2:ny]-demMask[, 1:(ny-1)])/dy

#########################################################
#Assign slopes for the upper and right border cells (i.e. where slopes can't be calcualted)
#Look for any upper and left borders where the slope couldn't be calculated and 
#repate the i-1 or j-1 slope
slopex2=slopex1
slopey2=slopey1
  
### Right  
borderR=slopeR=matrix(0, nrow=nx, ncol=ny)
  
#Turn the NA's into zeros for all the right border cells
borderR[1:(nx-1),]=(mask[1:(nx-1),]+mask[2:nx,])* mask[1:(nx-1),]
borderR[nx,]=1  
borderR=borderR*bordi
Rlist=which(borderR==1)
slopex2[Rlist]=0

#Find the right border cells that also have a cell to their left in the mask
borderR2=matrix(0, nrow=nx, ncol=ny)
borderR2[nx,]=1
borderR2[2:(nx-1),]=(mask[2:(nx-1),]+mask[3:nx,]) * mask[1:(nx-2),]
borderR2=borderR2*bordi #getting rid of any edge not on the border
borderR2[borderR2>1]=0

#if there is a cell to the left to use then fill in the slope for the border cell with this value
slopeR[2:nx,]=slopex2[1:(nx-1),]*borderR2[2:nx,]
slopeR[is.na(slopeR)==T]=0
slopex2=slopex2+slopeR
  
### Top
borderT=slopeT=matrix(0, nrow=nx, ncol=ny)
  
#Turn the NA's into zeros for all the right border cells
borderT[,1:(ny-1)]=(mask[,1:(ny-1)]+mask[,2:ny])* mask[,1:(ny-1)]
borderT[,ny]=1
borderT=borderT*bordi
Tlist=which(borderT==1)
slopey2[Tlist]=0
  
#Find the right border cells that also have a cell to their left in the mask
borderT2=matrix(0, nrow=nx, ncol=ny)
borderT2[,ny]=1
borderT2[,2:(ny-1)]=(mask[,2:(ny-1)]+mask[,3:ny]) * mask[,1:(ny-2)]
borderT2=borderT2*bordi #getting rid of any edge not on the border
borderT2[borderT2>1]=0
  
#if there is a cell to the left to use then fill in the slope for the border cell with this value
slopeT[,2:ny]=slopey2[,1:(ny-1)]*borderT2[,2:ny]
slopeT[is.na(slopeT)==T]=0
slopey2=slopey2+slopeT

######################################
#Make masks of the primary flow directions
#Make lists of the primary directions now that directions are all filled in
downlist=which(direction==d4[1])
leftlist=which(direction==d4[2])
uplist=which(direction==d4[3])
rightlist=which(direction==d4[4])

downlist.arr=which(direction==d4[1], arr.ind=T)
leftlist.arr=which(direction==d4[2], arr.ind=T)
uplist.arr=which(direction==d4[3], arr.ind=T)
rightlist.arr=which(direction==d4[4], arr.ind=T)

#masks of which cells contain primary flow direciton slope calcs
#Because the slopes are face centered - a primary flow direction up or right will mean that the slope for
#that cell is a primary slope
#However a flow direction down or left indicates that the slope  of the i-1 (or j-1) cell is a primary slope
#Thus some cells may have primary slopes in x and y while some may have them in neither even though every 
#grid cell only has one direction. 
ymask=xmask=matrix(0, nrow=nx, ncol=ny)
#mask of cells with primary flow in x and y direction, 
# signs indicate sign of the slope consistent - i.e. for flow in the positive x direciton (right) you need a negative slope
ymask[uplist]=-1
xmask[rightlist]=-1
  if(length(leftlist.arr>0)){
    for(ii in 1:nrow(leftlist.arr)){
	    xindex=max((leftlist.arr[ii,1]-1),1) #max statment means that if you have a left facing cell on the left border of the domain keep the slope at [i] as primary
	    if(mask[xindex, leftlist.arr[ii,2]]==0){xindex=xindex+1} #if the left cell falls outside the mask use the current cell for this border
	     xmask[xindex, leftlist.arr[ii,2]]=1
    }
  }
  if(length(downlist.arr>0)){
	  for(ii in 1:nrow(downlist.arr)){
	    yindex=max((downlist.arr[ii,2]-1),1) #max statment means that if you have a down facing cell on the lower border of the domain keep the slope at [i] as primary
	    if(mask[downlist.arr[ii,1],yindex]==0){yindex=yindex+1} #if the lower cell falls outside the mask use the current cell for this border
	    ymask[downlist.arr[ii,1], yindex]=1
	  }
  }
	ylist=which(ymask!=0) #primary flow direction y slope calcs
	xlist=which(xmask!=0) #primary flow direction x slope calcs

###################################	
#Do a check to see that slope direcitons are consistent with flow direcitons for primary flows
	fixPx=which(sign(slopex2)==-1 & xmask==1)
	    slopex2[fixPx]=abs(slopex2[fixPx])
	fixNx=which(sign(slopex2)==1 & xmask==-1)
	    slopex2[fixNx]=-abs(slopex2[fixNx])
	fixPy=which(sign(slopey2)==-1 & ymask==1)
	    slopey2[fixPy]=abs(slopey2[fixPy])
	fixNy=which(sign(slopey2)==1 & ymask==-1)
	    slopey2[fixNy]=-abs(slopey2[fixNy])
	Sinklist=c(fixPx, fixNx, fixPy, fixNy)
	
###################################
#If an lower limit on slopes is set (i.e. minslope is greater than zero)
# Then  apply the minimum slope threshold to any primary flow direction slope
#  while maintaining the direction consisten with the flow direction file
if(minslope>=0){
  print(paste("Limiting slopes to minimum", minslope))
  #x slopes
  xclipP=which(abs(slopex2)<minslope & xmask==1)	
  slopex2[xclipP]=minslope
	xclipN=which(abs(slopex2)<minslope & xmask==(-1))
	slopex2[xclipN]=(-minslope)
	  
	#y slopes
	yclipP=which(abs(slopey2)<minslope & ymask==1)
	slopey2[yclipP]=minslope
	yclipN=which(abs(slopey2)<minslope & ymask==(-1))
	slopey2[yclipN]=(-minslope)
	#print(paste('min adjustment', length(xclipP), range(slopex2[xclipP]), length(xclipN), length(yclipP), length(yclipN), range(slopey2[yclipP])))
}
	
###################################
#If an upper limit on slopes is set (i.e. maxslope is positive)
if(maxslope>=0){
	print(paste("Limiting slopes to maximum absolute value of", maxslope))
	#x slopes
	xclipP=which(slopex2>maxslope)
	slopex2[xclipP]=maxslope
	xclipN=which(slopex2<(-maxslope))
	slopex2[xclipN]=(-maxslope)

	#y slopes
	yclipP=which(slopey2>maxslope)
	slopey2[yclipP]=maxslope
	yclipN=which(slopey2<(-maxslope))
	slopey2[yclipN]=(-maxslope)
}


###################################
#If a maximum secondary/primary slope ratio is set (i.e. scalesecond==T)
if(secondaryTH>=0){
   if(secondaryTH==0){
     print(paste("Limiting the ratio of secondary to primary slopes", secondaryTH))
      slopex2[-xlist]=0
      slopey2[-ylist]=0
   } else{print("Options for nonzero secondary scaling not currently availabie please set  secondaryTH to -1 or 0  ")}
}



###################################
#Check for flat cells
# a cell is only flat if all 4 faces of the cell are flat
#flattest=matrix(0, ncol=ny, nrow=nx)
flattest=abs(slopex2[2:nx,2:ny])+abs(slopex2[1:(nx-1),2:ny])+abs(slopey2[2:nx,2:ny])+abs(slopey2[2:nx,1:(ny-1)])
nflat=length(which(flattest==0 & mask[2:nx, 2:ny]))
flats=which(flattest==0, arr.ind=T)+1
if(nflat>0){
  print(paste("WARNING:", nflat, "Flat cells found"))
}


###########################
#replace the NA's with 0s
nax=which(is.na(slopex2==T))
nay=which(is.na(slopey2==T))
slopex2[nax]=0
slopey2[nay]=0

output_list=list("slopex"=slopex2, "slopey"=slopey2, "direction"=direction, "Sinks"=Sinklist)
return(output_list)

} # end function
