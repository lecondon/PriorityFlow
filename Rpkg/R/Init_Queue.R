#' Initilze queue for topographic processing
#' 
#' Sets up a queue and initilizes marked and step matrices for DEM processing
#'
#' @inheritParams D4TraverseB
#' @param initmask Mask of the same dimensions of dem denoting a subset of cells to be considered for the queue (e.g. if you want to setup a run starting with only river cells). 
#' Note: if no init mask is included every border cell will be added to the queue
#' @param domainmask Mask of the domain extent to be considered. If no domain mask is provided boundaries will be calcualted from the rectangular extent
#' @param border Alternatively you can input your own border rathern than having it be calculated from the domain mask. For example if you want to have the river network and the borders combined you can input this as a border.
#' @export
InitQueue=function(dem, initmask, domainmask, border, d4=c(1,2,3,4)){
	
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

#initilize queue and matrices
ny=ncol(dem)
nx=nrow(dem)
queue=NULL
marked=matrix(0, nrow=nx, ncol=ny)
step=matrix(0, nrow=nx, ncol=ny)

#setup flow directions
#D4 neighbors
kd=matrix(0, nrow=4, ncol=3) #ordered down, left top right
kd[,1]=c(0,-1,0,1)
kd[,2]=c(-1,0,1,0)
kd[,3]=c(d4[1], d4[2], d4[3], d4[4])  


if(missing(initmask)){
	print("No init mask provided all border cells will be added to queue")
	initmask=matrix(1, nrow=nx, ncol=ny)
}

if(missing(domainmask)){
	print("No domain mask provided using entire domain")
	domainmask=matrix(1, nrow=nx, ncol=ny)
}


#Setup the border
if(missing(border)){
	print("No border provided, setting border using domain mask")
	border=matrix(1, nrow=nx, ncol=ny)
	border[2:(nx-1), 2:(ny-1)]= domainmask[1:(nx-2), 2:(ny-1)] + 													
	            domainmask[3:nx, 2:(ny-1)] +
							domainmask[2:(nx-1), 1:(ny-2)] +
							domainmask[2:(nx-1), 3:ny]
	border=border*domainmask
	border[which(border<4 & border!=0)]=1
	border[border==4]=0
}



basin=matrix(0, nrow=nx, ncol=ny)
maskbound=initmask*border
blist=which(maskbound>0) #array indices
binlist=which(maskbound>0, arr.ind=T) #xy indices
queue=cbind(binlist, dem[blist])
marked[blist]=1
basin[blist]=1:length(blist)

#assign flow direction to point out of the domain
direction=matrix(NA, nrow=nx, ncol=ny)
for(i in 1:nrow(queue)){
  xtemp=queue[i,1]
  ytemp=queue[i,2]
  temp=rep(0,4)
  for(d in 1:4){
    temp[d]=domainmask[xtemp+kd[d,1], ytemp+kd[d,2]]
  }
  dtemp=which.min(temp)
  direction[xtemp,ytemp]=kd[dtemp,3]
  
  
}

output_list=list("mask"=initmask,"queue" = queue, "marked"=marked, "basins"=basin, "direction"=direction)
return(output_list)
}
