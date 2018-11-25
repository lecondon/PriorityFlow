FindOrphan=function(dem, mask, marked){
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

#function to look for unprocessed river cells that have d8 neigbors on the river network
#or on the boundary

#matrix of offsets in x and y for D8 neigbors
kd=matrix(0, ncol=2, nrow=8)
kd[,1]=c(0, -1, -1, -1, 0, 1,1,1)
kd[,2]=c(-1, -1, 0, 1, 1, 1, 0, -1)

queue=NULL

t1=proc.time()
# Look for unprocessed river cells that have a d8 neigbhor on the river network
mlist=which(mask==1 & marked==0)
missed=matrix(0,nx,ny)
missed[mlist]=1
#print(paste(sum(missed), "channel cells missed"))

#Check and see if missed cells have D8 neigbors
#that were marked in the last traverse
#count the number of marked cells by summing the marked mask shifted for each of the d8 directions
ncount=matrix(0,nx,ny)
ncount[2:(nx-1), 2:(ny-1)] = marked[3:nx, 3:ny] +
						 	marked[1:(nx-2), 1:(ny-2)] +
						 	marked[3:nx, 1:(ny-2)] +
						 	marked[1:(nx-2), 3:ny] +
							marked[2:(nx-1), 3:ny] +
							marked[2:(nx-1), 1:(ny-2)] +
						    marked[3:nx, 2:(ny-1)] +
						    marked[1:(nx-2), 2:(ny-1)]
ncount=ncount*missed
norphan=length(which(ncount>0)) #number of orphaned branches: i.e. missed grid cells where the number of d8 marked neigbors >0


if(norphan>0){
	#add the marked neigbors of any orphan cell to a new queue
	addloc=which(ncount>0, arr.ind=T) # get the xy locations of the  orphans

	for(n in 1:norphan){
		xn=addloc[n,1]
		yn=addloc[n,2]
		#for each orphan look in all 8 directions and add any marked cell to the queue
		for(k in 1:8){
			xtemp=xn+kd[k,1]
			ytemp=yn+kd[k,2]
			if(marked[xtemp,ytemp]==1){
				queue=rbind(queue, c(xtemp, ytemp,
				dem[xtemp,ytemp]))
			} #end if
		} #end for k in 1:8
	} #end for n in 1:norphan
} else{
	print("No Orphans Found")
}

output_list=list("norphan"=norphan,"queue" = queue)
return(output_list)

}
