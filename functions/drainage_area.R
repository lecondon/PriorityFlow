drainageArea=function(direction, mask, d4=c(1,2,3,4), printflag=F){
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

#DraingeArea - function to calculate numbe of cells draining to any cell
#given a flow direction file

#Mandatory Inputs:
# 1. direction: numerical matrix of d4 flow directions

#Optional Inputs:
# 1. d4: directional numbering system for the direction matrix provided
#	order must be: down, left, top,right
#   defaults to: 1,2,3,4
# 2. mask: Mask with ones for cells to be processed and zeros for everything else - defaults to a mask of all 1's

nx=nrow(direction)
ny=ncol(direction)

if(missing(mask)){mask=matrix(1, nrow=nx, ncol=ny)} #default to processing everything

#Setup the border
border=matrix(1, nrow=nx, ncol=ny)
border[2:(nx-1), 2:(ny-1)]= mask[1:(nx-2), 2:(ny-1)] + 													mask[3:nx, 2:(ny-1)] +
							mask[2:(nx-1), 1:(ny-2)] +
							mask[2:(nx-1), 3:ny]
border=border*mask
border[which(border<4 & border!=0)]=1
border[border==4]=0


#set the borders to be the edge of the domain
#if(missing(border)){
#	border=matrix(0, nrow=nx, ncol=ny)
#	border[1:nx,c(1,ny)]=1
#	border[c(1,nx), 1:ny]=1
#}


#initilize drinage area matrix
drainarea=matrix(1, nrow=nx, ncol=ny)

#D4 neighbors
kd=matrix(0, nrow=4, ncol=2) #ordered down, left top right
kd[,1]=c(0,-1,0,1)
kd[,2]=c(-1,0,1,0)

#make masks of which cells drain down, up, left right
down=up=left=right=matrix(0, nrow=nx, ncol=ny)
down[which(direction==d4[1])]=1
left[which(direction==d4[2])]=1
up[which(direction==d4[3])]=1
right[which(direction==d4[4])]=1

#calculate the number of cells draining to any cell
draincount=matrix(0, nrow=nx, ncol=ny)
draincount[,1:(ny-1)]=draincount[,1:(ny-1)]+down[,2:ny]
draincount[,2:ny]=draincount[,2:ny]+up[,1:(ny-1)]
draincount[1:(nx-1),]=draincount[1:(nx-1),]+left[2:nx,]
draincount[2:nx, ]=draincount[2:nx,]+right[1:(nx-1),]

#give values outside the mask and on the bordr a negative count so they aren't processed
draincount[which(mask==0)]=(-99)
#draincount[which(border==1)]=(-99) #commented this out becuase it messes up when border cells point in. Not sure if I need this for something else though

#initilize a queue with all the headwater cells (i.e. cells with zero upstream neigbors)
draintemp=draincount
queue=which(draintemp==0, arr.ind=T)
qlist=which(draintemp==0)
blist=cbind(which(draintemp>0), which(draintemp>0, arr.ind=T))
nqueue=nrow(queue)
#image(draincount)
#image(draincount, zlim=c(0,0.5))

ii=1

while(nqueue>0){
	if(printflag){print(paste("lap", ii, "ncell", nqueue))}

	#loop through the queue
	for(i in 1:nqueue){
		#look downstream add 1 to the area and subtract 1 from the drainage #
		xtemp=queue[i,1]
		ytemp=queue[i,2]

		#if its has a flow direction
		if(is.na(direction[xtemp,ytemp])==F){
			dirtemp=which(d4==direction[xtemp,ytemp])
			xds=xtemp+kd[dirtemp,1]
			yds=ytemp+kd[dirtemp,2]

			#add one to the area of the downstream cell as long as that cell is in the domain
			if(xds<=nx & xds>=1 & yds<=ny & yds>=1){
				drainarea[xds, yds]=drainarea[xds, yds]+drainarea[xtemp,ytemp]

			#subtract one from the number of upstream cells from the downstream cell
			draintemp[xds,yds]= draintemp[xds,yds] - 1
		} #end if in the domain extent
	} #end if not na

		#set the drain temp to -99 for current cell to indicate its been done
		draintemp[xtemp,ytemp]=-99
	} #end for i in 1:nqueue

	#make a new queue with the cells with zero upstream drains left
	#queue=which(draintemp==0, arr.ind=T)
	ilist=which(draintemp[blist[,1]]==0)
	queue=blist[ilist,2:3]
	if(length(ilist)!=length(blist)/3){
		blist=blist[-ilist,]
	} else{
		blist=NULL
		if(printflag){print("blist empty")}
	}

	nqueue=length(queue)/2

	if(nqueue==1){queue=matrix(queue, ncol=2, nrow=1)}
	if(length(blist)/3==1){blist=matrix(blist, ncol=3, nrow=1)}
	ii=ii+1
	#print(length(which(draintemp<0 & draintemp>(-99))))
}

drainarea=drainarea*mask
return(drainarea)
#par(mfrow=c(1,2))
#image.plot(drainarea)
#image.plot(draintemp, zlim=c(-1,4))
#length(which(draintemp==0)

}
