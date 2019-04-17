SlopeCalcUP=function(dem, direction, dx, dy,  mask, borders, borderdir=1,  d4=c(1,2,3,4),  minslope=1e-5, maxslope=-1, secondaryTH=-1, river_method=0, rivermask, subbasins, printflag=F, upflag=T){
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

#SlopeCalcUP - Function to calculate the slopes from a DEM
#Mandatory Inputs:
# dem - matrix of elevation values
# direction - flow directions from processing
# dx, dy - grid size


# Optional Inputs:
# mask - matrix with 1 for cells within the domain and 0 for cells outside the domain
# borders - a matrix with 1's for borders cells to default to pointing out, 2 to default to pointing in, and 0 for all non-border cells.

#borderddir - default value for border cells 1=point out, 2= point in

# d4: directional numbering system: the numbers you want to assign to down, left, top,right (defaults to 1,2,3,4)

# minslope: Minimum absolute slope value to apply to flat cells if needed. Defaults to 1e-5
# maxslope - Maximum absolute value of slopes. If this is set to -1 the slopes will not be limited. Default value is -1

# secondary threshold - maximum ratio of |secondary|/|primary| to be enforced.
# NOTE - this scaling occurs after any max threholds are applied. If this is set to -1 no scaling will be applied.

#rivermethod- Optional method to treat river cells differently from the rest of the domain
	# 0: default value, no special treatment for river cells
	# 1: Scale secondary slopes to zero along the river (requries river mask)
	# 2: Apply watershed mean slope to each river reach (requires river mask and subbasins)
	# 3: Apply the stream mean slope to each reach (requires river mask and subbasins)
	# NOTE: the river mask can be different from the rivers that were used to create the subbasins if desired (i.e. if you want to use a threshold of 100 to create subbasins but then apply to river cells with a threshold of 50)

#rivermask - Mask with 1 for river cells and 0 for other cells

#upflag - A flag indicating whether slope calc should be upwinded, defaults to T generating slopes that are consistent with ParFlow. If set to F
#    then all slopes will be calcualted as [i+1]-[i]


###############################################################
ny=ncol(dem)
nx=nrow(dem)

#If no mask is provided default to the rectugular domain
if(missing(mask)){
	print("No mask provided, initializing mask for complete domain")
	mask=matrix(1, nrow=nx, ncol=ny)
	borders=matrix(0, nrow=nx, ncol=ny)
	borders[,c(1,ny)]=1*borderdir
	borders[c(1,nx),]=1*borderdir
}

#if no border is provided create a border with everything pointing in or out according to borderdir
if(missing(borders)){
	borders=matrix(1, nrow=nx, ncol=ny)
	borders[2:(nx-1), 2:(ny-1)]= mask[1:(nx-2), 2:(ny-1)] + 																mask[3:nx, 2:(ny-1)] +
								mask[2:(nx-1), 1:(ny-2)] +
								mask[2:(nx-1), 3:ny]
	borders=borders*mask
	borders[which(borders<4 & borders!=0)]=1 *borderdir
	borders[borders==4]=0
}

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
slopex1[nx,]=slopex1[(nx-1),]
slopey1[,1:(ny-1)]=(demMask[,2:ny]-demMask[, 1:(ny-1)])/dy
slopey1[,ny]=slopey1[,(ny-1)]

#Assign slopes based on upwinding for all non border cells
if(upflag==T){
	print("upwinding slopes")
slopex2=slopex1
slopey2=slopey1
for(j in 2:(ny-1)){
	#print(j)
	for(i in 2:(nx-1)){
		#print(i)
		if(mask[i,j]==1 & borders[i,j]==0){

		# X-Direction
		###################
		#SINK: If slopex[i-1]<0 & slopex[i]>0
		#Assign slope=0
		if(slopex1[(i-1),j]<0 & slopex1[i,j]>0){
			slopex2[i,j]=0
			if(direction[i,j]==2 | direction[i,j]==4){
				print(paste("Problem! Local X sink found in primary direction. x=", i, " j=", j, sep=""))
			}
		}

		#Local Maximum: If slopex[i-1]>0 & slopex[i]<0
		# If the primary flow direction is in the x
		# Then assign the slope consistent with this.
		# If not, choose the maximum absolute value
		if(slopex1[(i-1),j]>=0 & slopex1[i,j]<=0){
			#if the primary direction is left use the i-1 slope
			if(direction[i,j]==d4[2]){
				slopex2[i,j]=slopex1[(i-1),j]
			} else if(direction[i,j]==d4[4]){
			#if the primary direction is right use the i slope
				slopex2[i,j]=slopex1[i,j]
			} else{
				#If the primary direction is in the y then just choose the steepest one
				if(abs(slopex1[(i-1),j])> abs(slopex1[i,j])){
					slopex2[i,j]=slopex1[(i-1),j]
				}
			}
		}

		#Flow Through from right: If slopex[i-1]>0 & slopex[i]>0
		if(slopex1[(i-1),j]>=0 & slopex1[i,j]>=0){
			slopex2[i,j]=slopex1[(i-1),j]
		}

		#NOTE- Nothing to do for flow through from left
		#this would get assigned slopex[i,j] which is how slopex2 is initialized

		# Y-Direction
		###################
		#SINK: If slopey[j-1]<0 & slopey[j]>0
		#Assign slope=0
		if(slopey1[i,(j-1)]<0 & slopey1[i,j]>0){
			slopey2[i,j]=0
			if(direction[i,j]==1 | direction[i,j]==3){
				print(paste("Problem! Local Y sink found in primary direction. x=", i, " j=", j, sep=""))
			}
		}

		#Local Maximum: If slopey[j-1]>0 & slopex[j]<0
		# If the primary flow direction is in the y
		# Then assign the slope consistent with this.
		# If not, choose the maximum absolute value
		if(slopey1[i,(j-1)]>=0 & slopey1[i,j]<=0){
			#if the primary direction is down use the j-1 slope
			if(direction[i,j]==d4[1]){
				slopey2[i,j]=slopey1[i,(j-1)]
			} else if(direction[i,j]==d4[3]){
			#if the primary direction is up use the i slope
				slopey2[i,j]=slopey1[i,j]
			} else{
				#If the primary direction is in the x then just choose the steepest one
				if(abs(slopey1[i,(j-1)])> abs(slopey1[i,j])){
					slopey2[i,j]=slopex1[i,(j-1)]
				}
			}
		}

		#Flow Through from top: If slopex[j-1]>0 & slopex[j]>0
		if(slopey1[i,(j-1)]>=0 & slopey1[i,j]>=0){
			slopey2[i,j]=slopey1[i,(j-1)]
		}

		#NOTE- Nothing to do for flow through from bottom
		#this would get assigned slopey[i,j] which is how slopey2 is initialized

	}#end if inside the mask and not a border cell
	} #end for i
} #end for j
 #end if upflag
} else{
	#If you are not using the upwinded slopes just use the [i+1]-i Calculations
  print("standard slope calc")
	slopex2=slopex1
	slopey2=slopey1
}
###################
#Assign flow directions and slopes for border cells
bordi=borders #making a border indicator file with all 1's for any border cell
bordi[bordi>1]=1
interior=mask-bordi #mask of non boundary cells inside the domain

#Find the left border cells
borderL=matrix(0, nrow=nx, ncol=ny)
borderL[1,]=1
borderL2=borderL
#looking for cells where mask[i-1]=0 and mask[i]=1 (mask[2:nx]+mask[1:nx-1]=1)
#filtering out cells where there isn't a cell to the right to calculat the slope from (i.e. *mask[3:nx])
borderL[2:(nx-1),]=(mask[2:(nx-1),]+mask[1:(nx-2),]) * interior[3:nx,]
borderL=borderL*bordi #getting rid of any edge not on the border
borderL[borderL>1]=0
#second verions allows for any neighbor in the maks for slope calc (i.e. not just non border cells)
borderL2[2:(nx-1),]=(mask[2:(nx-1),]+mask[1:(nx-2),]) * mask[3:(nx),]
borderL2=borderL2*bordi #getting rid of any edge not on the border
borderL2[borderL2>1]=0
#left borders get value of 2
borderL=borderL*2
borderL2=borderL2*2

#Find the right border cells
borderR=matrix(0, nrow=nx, ncol=ny)
borderR[nx,]=1
borderR2=borderR
#looking for cells where mask[i+1]=0 and mask[i]=1 (mask[2:nx]+mask[1:nx-1]=1)
#filtering out cells where there isn't a cell to the left to calculat the slope from (i.e. *mask[1:(nx-2),])
borderR[2:(nx-1),]=(mask[2:(nx-1),]+mask[3:nx,]) * interior[1:(nx-2),]
borderR=borderR*bordi #getting rid of any edge not on the border
borderR[borderR>1]=0
#second verions allows for any neighbor in the maks for slope calc (i.e. not just non border cells)
borderR2[2:(nx-1),]=(mask[2:(nx-1),]+mask[3:nx,]) * mask[1:(nx-2),]
borderR2=borderR2*bordi #getting rid of any edge not on the border
borderR2[borderR2>1]=0
#right borders get value of 4
borderR=borderR*4
borderR2=borderR2*4

#Find the lower (bottom) border cells
borderB=matrix(0, nrow=nx, ncol=ny)
borderB[,1]=1
borderB2=borderB
#looking for cells where mask[j-1]=0 and mask[j]=1 (mask[2:ny]+mask[1:ny-1]=1)
#filtering out cells where there isn't a cell above to calculat the slope from
borderB[,2:(ny-1)]=(mask[, 2:(ny-1)]+mask[,1:(ny-2)]) * interior[,3:ny]
borderB=borderB*bordi #getting rid of any edge not on the border
borderB[borderB>1]=0
#second verions allows for any neighbor in the maks for slope calc (i.e. not just non border cells)
borderB2[,2:(ny-1)]=(mask[, 2:(ny-1)]+mask[,1:(ny-2)]) * mask[,3:ny]
borderB2=borderB2*bordi #getting rid of any edge not on the border
borderB2[borderB2>1]=0
#bottom borders get value of 1
borderB=borderB*1
borderB2=borderB2*1

#Find the upper (top) border cells
borderT=matrix(0, nrow=nx, ncol=ny)
borderT[,ny]=1
borderT2=borderT
#looking for cells where mask[j+1]=0 and mask[j]=1 (mask[2:ny]+mask[1:ny-1]=1)
#filtering out cells where there isn't a cell below to calculat the slope from
borderT[,2:(ny-1)]=(mask[,2:(ny-1)]+mask[,3:ny]) * interior[, 1:(ny-2)]
borderT=borderT*bordi #getting rid of any edge not on the border
borderT[borderT>1]=0
#second verions allows for any neighbor in the mask for slope calc (i.e. not just non border cells)
borderT2[,2:(ny-1)]=(mask[,2:(ny-1)]+mask[,3:ny]) * mask[, 1:(ny-2)]
borderT2=borderT2*bordi #getting rid of any edge not on the border
borderT2[borderT2>1]=0
#top borders get value of 3
borderT=borderT*3
borderT2=borderT2*3

#assign the edge face to be calculate as the max of the borders. In most cases there is just one
#non zero border so this is not really a choice
#In cases where there is more than one boundary it will have an arbitrary preference in the following order - Right, Top, Left, Bottom
borlist=which(bordi==1, arr.ind=T)
borlisti=which(bordi==1)
bordsum=(borderB+borderT+borderL+borderR)
borddir=matrix(0, nrow=nx, ncol=ny)
for(k in 1:nrow(borlist)){
	i=borlist[k,1]
	j=borlist[k,2]
		#If you can choose a border edge with an internal cell opposite do this first
		if(bordsum[i,j]>0){
			borddir[i,j]=max(borderB[i,j],borderT[i,j],borderL[i,j],borderR[i,j])
		#If not resort to picking an edge with a boundary cell opposite (this mostly happens on corners where both the d4 neigbors are also borders)
		} else{
			borddir[i,j]=max(borderB2[i,j],borderT2[i,j],borderL2[i,j],borderR2[i,j])
		}
}

test=length(which((borddir-bordi)<0))
missinglist=which((borddir-bordi)<0)
missinglistA=which((borddir-bordi)<0, arr.ind=T)
if(test>0){print(paste("ERROR:", test, "border cells missed!!"))}
#these errors occur with sinks or some special cases where there arent ANY
#usable neighbors e.g. the center cell of the example
#[,1] [,2] [,3] [,4] [,5]
#[1,]    1    0    0    0    0
#[2,]    1    0    0    0    0
#[3,]    1    1    1    1    0
#[4,]    1    1    0    0    0
#[5,]    1    0    0    0    0
#These cells get assigned an arbitrary direction and minimum slope below


#calculate the number of cells draining to any cell
down=up=left=right=matrix(0, nrow=nx, ncol=ny)
down[which(direction==d4[1])]=1
left[which(direction==d4[2])]=1
up[which(direction==d4[3])]=1
right[which(direction==d4[4])]=1
draincount=matrix(0, nrow=nx, ncol=ny)
draincount[,1:(ny-1)]=draincount[,1:(ny-1)]+down[,2:ny]
draincount[,2:ny]=draincount[,2:ny]+up[,1:(ny-1)]
draincount[1:(nx-1),]=draincount[1:(nx-1),]+left[2:nx,]
draincount[2:nx, ]=draincount[2:nx,]+right[1:(nx-1),]

#if a cell on the boder has a another interior cell draining to it then it should be pointing out.
drainout=which(draincount[borlisti]>0)
borders[borlisti[drainout]]=1

#Calculate the slopes
for(k in 1:nrow(borlist)){
	i=borlist[k,1]
	j=borlist[k,2]

	#border on the top
	if(borddir[i,j]==3){
		slopex2[i,j]=0
		#if its pointing out make slopey negative
		if(borders[i,j]==1){
			slopey2[i,j]=-abs(slopey1[i,(j-1)])
			direction[i,j]=3
		}else{
			#if its pointing in make slopey positive
			slopey2[i,j]=abs(slopey1[i,(j-1)])
			direction[i,j]=1
		}
	}

	#border on the right
	if(borddir[i,j]==4){
		slopey2[i,j]=0
		#if its pointing out make slopex negative
		if(borders[i,j]==1){
			slopex2[i,j]=-abs(slopex1[(i-1),j])
			direction[i,j]=4
		}else{
			#if its pointing in make slopex positive
			slopex2[i,j]=abs(slopex1[(i-1),j])
			direction[i,j]=2
		}
	}

	#border on the bottom
	if(borddir[i,j]==1){
		slopex2[i,j]=0
		#if its pointing out make slopey positive
		if(borders[i,j]==1){
			slopey2[i,j]=abs(slopey1[i,j])
			direction[i,j]=1
		}else{
			#if its pointing in make slopey negative positive
			slopey2[i,j]=-abs(slopey1[i,j])
			direction[i,j]=3
		}
	}

	#border on the left
	if(borddir[i,j]==2){
		slopey2[i,j]=0
		#if its pointing out make slopex positive
		if(borders[i,j]==1){
			slopex2[i,j]=abs(slopex1[i,j])
			direction[i,j]=2
		}else{
			#if its pointing in make slopex negative
			slopex2[i,j]=-abs(slopex1[i,j])
			direction[i,j]=4
		}
	}
}

#filling in the missing values with an arbitrary direciotn
#need to circle back to this
direction[missinglist]=4
slopex2[missinglist]=-minslope
slopey2[missinglist]=0

#Make lists of the primary directions now that directions are all filled in
downlist=which(direction==d4[1])
leftlist=which(direction==d4[2])
uplist=which(direction==d4[3])
rightlist=which(direction==d4[4])
ylist=c(uplist, downlist) #list of cells with primary flow in the x direction
xlist=c(rightlist, leftlist) #list of cells with primary flow in the y direction
ymask=xmask=matrix(0, nrow=nx, ncol=ny)#mask of cells with primary flow in x and y direction, signs indicate direction
ymask[downlist]=1
ymask[uplist]=-1
xmask[leftlist]=1
xmask[rightlist]=-1

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
	print(paste("Limiting the ratio of secondary to primary slopes", secondaryTH))

	#Make matrices of primary and secondary slopes and ratios
	primary=abs(slopex2)
	primary[ylist]=abs(slopey2[ylist])
	secondary=abs(slopey2)
	secondary[ylist]=abs(slopex2[ylist])
	ratio=secondary/primary

	scalelist=which(ratio>secondaryTH)
	for(i in 1:length(scalelist)){
		temp=scalelist[i]
		#if primary direction is in x = then scale the y slopes
		if(direction[temp]==d4[2] || direction[temp]==d4[4]){
			slopey2[temp]=sign(slopey2[temp])*secondaryTH*abs(slopex2[temp])
		}
		#if primary direction is in y = then scale the x slopes
		if(direction[temp]==d4[1] || direction[temp]==d4[3]){
			slopex2[temp]=sign(slopex2[temp])*secondaryTH*abs(slopey2[temp])
		}
	}
}

###################################
#Separate processing for the river cells
#Option 1: Just turn off secondary slopes in the river cells
if(river_method==1){
	if(missing(rivermask)){
		print("WARNING: No rivermask provided, slopes not adjusted")
	}else{
		print("River Method 1: setting secondary slopes to zero along the river")
		print("Removing secondary slopes along river mask")
		#zero out the slopes in the secondary directions over the whole domain
		xtemp=slopex2
		xtemp[ylist]=0
		ytemp=slopey2
		ytemp[xlist]=0

		#merge in these slopes over the river mask
		rivlist=which(rivermask==1)
		slopex2[rivlist]=xtemp[rivlist]
		slopey2[rivlist]=ytemp[rivlist]
}
}

# Option 2: Assign the subbasin mean slope to the river cells in the primary direciton and set the secondary direction to 0
if(river_method==2){
	print("River Method 2: assigning average watershed slope to river cells by watershed")
	nbasin=max(subbasins)
	savg=rep(0,nbasin)
	count=rep(0,nbasin)
	#calculate the subbasin average slope
	#NEED to sort out for x list and ylist so just averaging the primary
	for(i in 1:nx){
		for(j in 1:ny){
			if(subbasins[i,j]>0){
				savg[subbasins[i,j]]=savg[subbasins[i,j]]+abs(slopex2[i,j])*abs(xmask[i,j])+abs(slopey2[i,j])*abs(ymask[i,j])
				count[subbasins[i,j]]=count[subbasins[i,j]] + 1
			} #end if
		} #end for j
	}#end for i
	savg=savg/count

	rivlist=which(rivermask==1)
	nriv=length(rivlist)
	#fill in the river slopes with the average for their subbasin
	for(i in 1:nriv){
		rtemp=rivlist[i]
		sbtemp=subbasins[rtemp]
		#print(paste(i, rtemp, sbtemp))
		if(sbtemp>0){
			slopex2[rtemp]=savg[sbtemp]*xmask[rtemp]
			slopey2[rtemp]=savg[sbtemp]*ymask[rtemp]
		} #end if
	}#end for i in 1:nriv
} # end if river method ==2

# Option 3: Assign the subbasin mean river slope to the river cells in the primary direciton and set the secondary direction to 0
if(river_method==3){
	print("River Method 3: assigning average river slope to river cells by watershed")
	nbasin=max(subbasins)
	savg=rep(0,nbasin)
	count=rep(0,nbasin)
	#calculate the subbasin average slope
	#NEED to sort out for x list and ylist so just averaging the primary
	for(i in 1:nx){
		for(j in 1:ny){
			if(subbasins[i,j]>0 & rivermask[i,j]==1){
				savg[subbasins[i,j]]=savg[subbasins[i,j]]+abs(slopex2[i,j])*abs(xmask[i,j])+abs(slopey2[i,j])*abs(ymask[i,j])
				count[subbasins[i,j]]=count[subbasins[i,j]] + 1
			} #end if
		} #end for j
	}#end for i
	savg=savg/count

	rivlist=which(rivermask==1)
	nriv=length(rivlist)
	#fill in the river slopes with the average for their subbasin
	for(i in 1:nriv){
		rtemp=rivlist[i]
		sbtemp=subbasins[rtemp]
		#print(paste(i, rtemp, sbtemp))
		if(sbtemp>0){
			slopex2[rtemp]=savg[sbtemp]*xmask[rtemp]
			slopey2[rtemp]=savg[sbtemp]*ymask[rtemp]
		} #end if
	}#end for i in 1:nriv
} # end if river method ==3

###################################
#Check for flat cells
nflat=length(which(slopex2==0 & slopey2==0))
if(nflat!=0){
	print(paste("WARNING:", nflat, "Flat cells found"))
	flatloc=which(slopex2==0 & slopey2==0, arr.ind=T)
	flatlist=which(slopex2==0 & slopey2==0)
	if(printflag){
		print("Flat locations (note this is x,y)")
		print(flatloc)
	}
	#impose a minimum slope in the primary direction for flat cells
	for(i in 1:nflat){
		dtemp=direction[flatlist[i]]
		if(dtemp==d4[1]){slopey2[flatlist[i]]=minslope} #down
		if(dtemp==d4[2]){slopex2[flatlist[i]]=minslope} #left
		if(dtemp==d4[3]){slopey2[flatlist[i]]=-minslope} #up
		if(dtemp==d4[4]){slopex2[flatlist[i]]=-minslope} #right
	}
	nflat=length(which(slopex2==0 & slopey2==0))
	print(paste("After processing:", nflat, "Flat cells left"))
}

###################################
#If an lower limit on slopes is set (i.e. minslope is greater than zero)
if(minslope>=0){
	print(paste("Limiting slopes to minimum", minslope))
	#x slopes
	xclipP=which(slopex2<minslope & slopex2>0 & xmask==1)
	slopex2[xclipP]=minslope
	xclipN=which(slopex2>(-minslope) & slopex2<0 & xmask==1)
	slopex2[xclipN]=(-minslope)

	#y slopes
	yclipP=which(slopey2<minslope & slopey2>0 & ymask==1)
	slopey2[yclipP]=minslope
	yclipN=which(slopey2>(-minslope) & slopey2<0 & ymask==1)
	slopey2[yclipN]=(-minslope)
	#print(paste('min adjustment', length(xclipP), range(slopex2[xclipP]), length(xclipN), length(yclipP), length(yclipN), range(slopey2[yclipP])))
}


###########################
#replace the NA's with 0s
nax=which(is.na(slopex2==T))
nay=which(is.na(slopey2==T))
slopex2[nax]=0
slopey2[nay]=0

output_list=list("slopex"=slopex2, "slopey"=slopey2, "direction"=direction)
return(output_list)

} # end function
