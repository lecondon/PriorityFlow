#'Calculate the slopes from a DEM
#'
#'This function will calculate slopes using standard or upwinding options and apply a range of smoothing options

#' @inheritParams D4TraverseB
#' @param  dx,dy Lateral grid cell resolution
#' @param borders A matrix with 1's for borders cells to default to pointing out, 2 to default to pointing in, and 0 for all non-border cells.
#' @param borderdir Default value for border cells 1=point out, 2= point in
#' @param adjust_borders Flag determining whether the border cells should be ajusted or not to point in a  specific direction
#' @param minslope Minimum absolute slope value to apply to flat cells if needed. Defaults to 1e-5
#' @param subbasins An optional matrix of subbasin values
#' @param maxslope  Maximum absolute value of slopes. If this is set to -1 the slopes will not be limited. Default value is -1
#' @param secondaryTH  secondary threshold - maximum ratio of |secondary|/|primary| to be enforced. NOTE - this scaling occurs after any max threholds are applied. If this is set to -1 no scaling will be applied.
#' @param river_secondaryTH  secondary threshold  to apply to the river cells if river method 1-3 is chosen- maximum ratio of |secondary|/|primary| to be enforced. NOTE - this scaling occurs after any max threholds are applied and the river smoothing is done if you are usnig optons 2 or 3.  The dafault value is zero. NOTE- there is not currently a '-1' equivalend for this threshold. secondary slopes along the river cells must either be scaled to 0 or some ratio of the primary slopes if one of the river methods is chosen.
#' @param river_method Optional method to treat river cells differently from the rest of the domain
#' @param rivermask Mask with 1 for river cells and 0 for other cells
#' @param printflag Print function progress
#' @param upflag A flag indicating whether slope calc should be downwinded to be consitent with the upwinding in the ParFlow OverlandFlow Boundary condition it defaults to T. If set to F then all slopes will be calcualted as [i+1]-[i] with no adjustments. This approach is consistent with the OverlandKin boundary condition in  ParFlow.

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
SlopeCalStan=function(dem, direction, dx, dy,  mask, borders, borderdir=1, adjust_borders=F, d4=c(1,2,3,4),  minslope=0, maxslope=-1, secondaryTH=-1, river_method=0, river_secondaryTH=0, rivermask, subbasins, printflag=F, upflag=F){
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
	borders[,c(1,ny)]=1*borderdir
	borders[c(1,nx),]=1*borderdir
}

#if no border is provided create a border with everything pointing in or out according to borderdir
if(missing(borders)){
	borders=matrix(1, nrow=nx, ncol=ny)
	borders[2:(nx-1), 2:(ny-1)]= mask[1:(nx-2), 2:(ny-1)] + 																
	              mask[3:nx, 2:(ny-1)] +
								mask[2:(nx-1), 1:(ny-2)] +
								mask[2:(nx-1), 3:ny]
	borders=borders*mask
	borders[which(borders<4 & borders!=0)]=1 *borderdir
	borders[borders==4]=0
}

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
#slopex1[nx,]=slopex1[(nx-1),]
slopey1[,1:(ny-1)]=(demMask[,2:ny]-demMask[, 1:(ny-1)])/dy
#slopey1[,ny]=slopey1[,(ny-1)]

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
  for(ii in 1:nrow(leftlist.arr)){
	  xindex=max((leftlist.arr[ii,1]-1),1) #max statment means that if you have a left facing cell on the left border of the domain keep the slope at [i] as primary
	  if(mask[xindex, leftlist.arr[ii,2]]==0){xindex=xindex+1} #if the left cell falls outside the mask use the current cell for this border
	   xmask[xindex, leftlist.arr[ii,2]]=1
	}
	for(ii in 1:nrow(downlist.arr)){
	  yindex=max((downlist.arr[ii,2]-1),1) #max statment means that if you have a down facing cell on the lower border of the domain keep the slope at [i] as primary
	  if(mask[downlist.arr[ii,1],yindex]==0){yindex=yindex+1} #if the lower cell falls outside the mask use the current cell for this border
	  ymask[downlist.arr[ii,1], yindex]=1
	}
	ylist=which(ymask!=0) #primary flow direction y slope calcs
	xlist=which(xmask!=0) #primary flow direction x slope calcs

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
	print(paste("Limiting the ratio of secondary to primary slopes", secondaryTH))

  if(upflag==T){
    #if you are using the cell centered slopes then its guranteed that every cell is 
    #either primary in x or primary in Y and the opposite direction is secondary
    #and you can take ratios of primary and secondary slopes on the same cell
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
  }else{
    #If you are using the face centered slopes this is not the case 
    #Refer to the calcualtions of the xlists and ylists above
   if(secondaryTH==0){
      slopex2[-xlist]=0
      slopey2[-ylist]=0
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
		print(paste("Scaling secondary slopes along river mask to", river_secondaryTH, "* primary slope"))

		#old approach before river_secondaryTH flag was added
		##zero out the slopes in the secondary directions over the whole domain
		#xtemp=slopex2
		#xtemp[ylist]=0
		#ytemp=slopey2
		#ytemp[xlist]=0

		#otherwise scale by whatever river secondary threshold is provided
		xtemp=slopex2
		xtemp_scaled=slopey2*river_secondaryTH #xslopes=yslopes*scaler
		xtemp[ylist]=xtemp_scaled[ylist] #where primary flow is in the y direction fill in x slopes with teh scaled values

		#repeat in y direction
		ytemp=slopey2
		ytemp_scaled=slopex2*river_secondaryTH #yslopes=xslopes*scaler
		ytemp[xlist]=ytemp_scaled[xlist]

		#merge in these slopes over the river mask
		rivlist=which(rivermask==1)
		slopex2[rivlist]=xtemp[rivlist]
		slopey2[rivlist]=ytemp[rivlist]

}
}

# Option 2: Assign the subbasin mean slope to the river cells in the primary direction and set the secondary direction to 0
if(river_method==2){
	print("River Method 2: assigning average watershed slope to river cells by watershed")
	print(paste("Scaling secondary slopes along river mask to", river_secondaryTH, "* primary slope"))

	nbasin=max(subbasins)
	savg=rep(0,nbasin)
	count=rep(0,nbasin)

	#Get the signs of the slopes before averaging to use in the secondary slope scaling
	xsign=sign(slopex2)
	ysign=sign(slopey2)
	xsign[xsign==0]=1 # for zero slopes default to positive
	ysign[ysign==0]=1 #for zero slopes default to positive

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
	for(b in 1:nbasin){
		savg[b]=max(minslope,savg[b]) #make sure the average slope is >=savg
	}
	savg[is.na(savg)]=minslope

	rivlist=which(rivermask==1)
	nriv=length(rivlist)
	#fill in the river slopes with the average for their subbasin
	for(i in 1:nriv){
		rtemp=rivlist[i]
		sbtemp=subbasins[rtemp]
		#print(paste(i, rtemp, sbtemp))
		if(sbtemp>0){
			#setting the primary slopes along the river (xmask and y masks are masks of primary flow direciton so secondary will be set to zero)
			slopex2[rtemp]=savg[sbtemp]*xmask[rtemp]
			slopey2[rtemp]=savg[sbtemp]*ymask[rtemp]

			#setting the secondary slopes along the river
			# By taking the inverse of the mask and adding the scale slope to the inital slope
			# This means it will add zero if its the primary flow direcito and it will add the scalaed secondary if not
			slopex2[rtemp]=slopex2[rtemp] + savg[sbtemp]*abs(1-abs(xmask[rtemp]))*xsign[rtemp]*river_secondaryTH
			slopey2[rtemp]=slopey2[rtemp] + savg[sbtemp]*abs(1-abs(ymask[rtemp]))*ysign[rtemp]*river_secondaryTH

		} #end if
	}#end for i in 1:nriv
} # end if river method ==2

# Option 3: Assign the subbasin mean river slope to the river cells in the primary direciton and set the secondary direction to 0
if(river_method==3){
	print("River Method 3: assigning average river slope to river cells by watershed")
	print(paste("Scaling secondary slopes along river mask to", river_secondaryTH, "* primary slope"))

	nbasin=max(subbasins)
	savg=rep(0,nbasin)
	count=rep(0,nbasin)
	#Get the signs of the slopes before averaging to use in the secondary slope scaling
	xsign=sign(slopex2)
	ysign=sign(slopey2)
	xsign[xsign==0]=1 # for zero slopes default to positive
	ysign[ysign==0]=1 #for zero slopes default to positive

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
	for(b in 1:nbasin){
		savg[b]=max(minslope,savg[b]) #make sure the average slope is >=savg
	}
	savg[is.na(savg)]=minslope

	rivlist=which(rivermask==1)
	nriv=length(rivlist)
	#fill in the river slopes with the average for their subbasin
	for(i in 1:nriv){
		rtemp=rivlist[i]
		sbtemp=subbasins[rtemp]
		#print(paste(i, rtemp, sbtemp))
		if(sbtemp>0){
			#setting the primary slopes along the river (xmask and y masks are masks of primary flow direciton so secondary will be set to zero)
			slopex2[rtemp]=savg[sbtemp]*xmask[rtemp]
			slopey2[rtemp]=savg[sbtemp]*ymask[rtemp]

			#setting the secondary slopes along the river
			# By taking the inverse of the mask and adding the scale slope to the inital slope
			# This means it will add zero if its the primary flow direcito and it will add the scalaed secondary if not
			slopex2[rtemp]=slopex2[rtemp] + savg[sbtemp]*abs(1-abs(xmask[rtemp]))*xsign[rtemp]*river_secondaryTH
			slopey2[rtemp]=slopey2[rtemp] + savg[sbtemp]*abs(1-abs(ymask[rtemp]))*ysign[rtemp]*river_secondaryTH

		} #end if

	}#end for i in 1:nriv
} # end if river method ==3

###################################
#Check for flat cells
# Only look in the cell centered case
if(upflag==T){
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
} else{
  # a cell is only flat if all 4 faces of the cell are flat
  #flattest=matrix(0, ncol=ny, nrow=nx)
  flattest=abs(slopex2[2:nx,2:ny])+abs(slopex2[1:(nx-1),2:ny])+abs(slopey2[2:nx,2:ny])+abs(slopey2[2:nx,1:(ny-1)])
  nflat=length(which(flattest==0))
  flats=which(flattest==0, arr.ind=T)+1
  if(nflat>0){
    print(paste("WARNING:", nflat, "Flat cells found"))
  }
}

###################################
### MOVED THIS UP higher--- delete this later
#If an lower limit on slopes is set (i.e. minslope is greater than zero)
#if(minslope>=0){
#	print(paste("Limiting slopes to minimum", minslope))
#	#x slopes
#	xclipP=which(slopex2<minslope & slopex2>0 & xmask==1)
#	slopex2[xclipP]=minslope
#	xclipN=which(slopex2>(-minslope) & slopex2<0 & xmask==1)
#	slopex2[xclipN]=(-minslope)
##	#y slopes
#	yclipP=which(slopey2<minslope & slopey2>0 & ymask==1)
#	slopey2[yclipP]=minslope
#	yclipN=which(slopey2>(-minslope) & slopey2<0 & ymask==1)
#	slopey2[yclipN]=(-minslope)
#	#print(paste('min adjustment', length(xclipP), range(slopex2[xclipP]), length(xclipN), length(yclipP), length(yclipN), range(slopey2[yclipP])))
#}


###########################
#replace the NA's with 0s
nax=which(is.na(slopex2==T))
nay=which(is.na(slopey2==T))
slopex2[nax]=0
slopey2[nay]=0

output_list=list("slopex"=slopex2, "slopey"=slopey2, "direction"=direction)
return(output_list)

} # end function
