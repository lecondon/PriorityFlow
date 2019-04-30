DelinWatershed=function(outlets, direction, d4=c(1,2,3,4), printflag=F){
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

# DelinWatershed - function to define the watershed for a point or set of outlet points based on the flow direction file

#Mandatory Inputs:
# 1. x,y coordinated of the outlet points or points to mask upstream areas for if there is just one point this can be input as c(x,y), if there are multiple points, this should be a matrix with a separate row for each point
# 2. direction: numerical matrix of d4 flow directions

#Optional Inputs:
# 1. d4: directional numbering system for the direction matrix provided
#	order must be: down, left, top,right
#   defaults to: 1,2,3,4
# 2. printflag: Optional flag to print out the number of cells in the queue durring iterations, defaults to F

nx=nrow(direction)
ny=ncol(direction)

#initilize a matrix to store the mask
marked=matrix(0, nrow=nx, ncol=ny)

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

#intialized the queue with the outlet points
if(length(outlets)==2){
  queue=matrix(outlets, nrow=1,  ncol=2) #if there is only one point pair format it into a matrix
}else{
  queue=outlets
}

nqueue=nrow(queue)
count0=0
ii=1
while(nqueue>0){
	if(printflag){print(paste("lap", ii, "ncell", nqueue))}
	queue2=NULL

	#loop through the queue
	for(i in 1:nqueue){
		xtemp=queue[i,1]
		ytemp=queue[i,2]
		#add one to the subbasin area for the summary
		marked[xtemp,ytemp]=1


		#look for cells that drain to this cell
		for(d in 1:4){
			xus=xtemp-kd[d,1]
			yus=ytemp-kd[d,2]
			if(xus*yus>0 & xus<=nx & yus<=ny){
			  if(is.na(direction[xus,yus])==F & marked[xus,yus]==0){
					if(direction[xus,yus]==d4[d]){
					  #print(c(xus,yus))
						marked[xus,yus]=1 #add the upstream cell to the mask
						queue2=rbind(queue2, c(xus,yus)) # and then add the upstream cell to the queue
					} #end if pointing to cell
			  } #end if in the mask
			} # end if its in the domain bounds
		} #end direction loop
	}
	if(length(queue2)>=2){
		queue=queue2
		nqueue=nrow(queue)
		ii=ii+1
	} else{nqueue=0}
}

masklist=which(marked==1, arr.ind=T)
xrange=range(masklist[,1])
yrange=range(masklist[,2])

output_list=list("watershed"=marked, 'xrange'=xrange, 'yrange'=yrange)

return(output_list)

}
