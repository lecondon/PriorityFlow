CalcSubbasins=function(direction, area, mask, d4=c(1,2,3,4), riv_th=50, printflag=F, merge_th=0){
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

# CalcSubbasins - function to divide the domain into subbasins with individual stream segments

#Mandatory Inputs:
# 1. direction: numerical matrix of d4 flow directions
# 2. area: drainage areas for every cell

#Optional Inputs:
# 1. d4: directional numbering system for the direction matrix provided
#	order must be: down, left, top,right
#   defaults to: 1,2,3,4
# 2. mask: Mask with ones for cells to be processed and zeros for everything else - defaults to a mask of all 1e's
# 3. riv_th= threshold for the drainage area minimum used desigate cells as river cells, defaults to 50
# 4. merge_th=after all the subbasins have been defined subbasins with areas <merg_th will be combined with their downstream neighbors (Defaults to 0 which means no merging will take place)

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
subbasin=matrix(0, nrow=nx, ncol=ny)
marked=matrix(0, nrow=nx, ncol=ny)

#D4 neighbors
kd=matrix(0, nrow=4, ncol=2) #ordered down, left top right
kd[,1]=c(0,-1,0,1)
kd[,2]=c(-1,0,1,0)

#make a river mask based on the drainage area threshold
rivers=area
rivers[area<riv_th]=0
rivers[area>=riv_th]=1
#print(sum(rivers))
#image.plot(rivers)

#make masks of which cells drain down, up, left right
down=up=left=right=matrix(0, nrow=nx, ncol=ny)
down[which(direction==d4[1])]=1
left[which(direction==d4[2])]=1
up[which(direction==d4[3])]=1
right[which(direction==d4[4])]=1

#calculate the number of river cells draining to any cell
draincount=matrix(0, nrow=nx, ncol=ny)
draincount[,1:(ny-1)]=draincount[,1:(ny-1)]+down[,2:ny]*rivers[,2:ny]
draincount[,2:ny]=draincount[,2:ny]+up[,1:(ny-1)]*rivers[,1:(ny-1)]
draincount[1:(nx-1),]=draincount[1:(nx-1),]+left[2:nx,]*rivers[2:nx,]
draincount[2:nx, ]=draincount[2:nx,]+right[1:(nx-1),]*rivers[1:(nx-1),]

#Identify all the headwater cells
headwater=matrix(0, nrow=nx, ncol=ny)
headwater[which(draincount==0 & rivers==1)]=1

#image.plot(headwater*3+draincount, zlim=c(0.5,4))

#give values outside the mask and on the border a negative count so they aren't processed
marked[which(mask==0)]=1


#start with all the headwater cells (i.e. cells with zero upstream neigbors)
blist=cbind(which(headwater==1), which(headwater==1, arr.ind=T))
nheadwater=nrow(blist)
ends=rivers
ends[blist[,1]]=2
#marked[qlist]=1
#subbasin[blist[,1]]=1:nrow(blist) #initialized the subbasin numbers
#image.plot(headwater)

#Get just the river areas to use for this
rivarea=area*rivers

index=0
#subbasin[blist[1,1]]=index
#marked[blist[1,1]]=1

subbasin=matrix(0, nrow=nx, ncol=ny)
marked=matrix(0, nrow=nx, ncol=ny)
first=T
###1. walk down from every headwater marking stream segments
for(i in 1:nheadwater){
	xtemp=blist[i,2]
	ytemp=blist[i,3]
	active=T
	index=index+1
	subbasin[xtemp,ytemp]=index
	marked[xtemp,ytemp]=1
	#print(paste("Starting new Branch", i,  index, subbasin[xtemp,ytemp]))
	#image.plot(subbasin)
	summarytemp=c(index, xtemp, ytemp, rep(0,4))

	while(active==T){
		#get the direction and find downstream cell
		dirtemp=which(d4==direction[xtemp,ytemp])
		xds=xtemp+kd[dirtemp,1]
		yds=ytemp+kd[dirtemp,2]

		#if the downstream neigbor hasn't already been procesed and its in the domain
		if(xds*yds>0 & xds<=nx & yds<=ny){
		if(marked[xds,yds]==0 & mask[xds,yds]==1){
			#Check the area difference
			accum=area[xds,yds]-area[xtemp, ytemp]

			#if there is a tributary coming in then start a new segment
			if(accum>riv_th){
				#print('Here')
				summarytemp[4:5]=c(xtemp,ytemp)
				summarytemp[6]=index+1
				index=index+1
				ends[xtemp,ytemp]=3
				ends[xds,yds]=2
				if(first==T){
					summary=summarytemp
					first=F
				}else{
					summary=rbind(summary, summarytemp)
				}
				summarytemp=c(index, xtemp, ytemp, rep(0,4))

				#print(paste("new index", index))
			}

			#assign subbasin number to the downstream cell and mark it off
			#print("step")
			subbasin[xds,yds]=index
			marked[xds,yds]=1
			xtemp=xds
			ytemp=yds
			#print(paste(xtemp,ytemp, index))
		} else{
			#if the downstream neighbor has been processed then move on to the next headwater cell
			#print("Alredy Marked")
			active=FALSE
			ends[xtemp,ytemp]=3
			summarytemp[4:5]=c(xtemp,ytemp)
			summarytemp[6]=subbasin[xds,yds]
		} #end else
	 } else{
	 	#print("Outside the domain")
		active=FALSE
		ends[xtemp,ytemp]=3
		summarytemp[4:5]=c(xtemp,ytemp)
		summarytemp[6]=-1
	 }
	} # end while
	if(first==T){
		summary=summarytemp
		first=F
	}else{
		summary=rbind(summary, summarytemp)
	}
} # end for
rownames(summary)=NULL
colnames(summary)=c("Basin_ID","start_x","start_y", "end_x", "end_y", "Downstream_Basin_ID", "Area")

#image.plot(subbasin, zlim=c(17,18))
#image.plot(marked)

#test=marked-rivers
#range(test)


###2. Get the drainage basins for every segement
subbasinA=subbasin

#start a queue with all the cells in the river
queue=which(subbasin>0, arr.ind=T)
qlist=which(subbasin>0)
blist=cbind(which(subbasin>0), which(subbasin>0, arr.ind=T))

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
		sbtemp=subbasinA[xtemp,ytemp]
		summary[sbtemp,7]=summary[sbtemp,7]+1

		#look for cells that drain to this cell
		for(d in 1:4){
			xus=xtemp-kd[d,1]
			yus=ytemp-kd[d,2]
			if(xus*yus>0 & xus<=nx & yus<=ny){
				if(mask[xus,yus]==1 & subbasinA[xus,yus]==0){
					if(direction[xus,yus]==d4[d]){
						subbasinA[xus,yus]=subbasinA[xtemp,ytemp] #assign the subbasin number to the upstream cell
						queue2=rbind(queue2, c(xus,yus)) #add to the queue
					} #end if pointing to cell
				} #end if its on the mask
			} # end if its in the domain bounds
		} #end direction loop
	}
	if(length(queue2)>=2){
		queue=queue2
		nqueue=nrow(queue)
		ii=ii+1
	} else{nqueue=0}
}


###3. if merge_th >0 look for basins with areas less than the merge threshold, that don't drain out of the domain and merge with their downstream neighbors

delete=NULL #list of subbasins to delete from summary list
if(merge_th>0){
	nsb=nrow(summary)
	for(i in 1:nsb){
		#check if area is less than the threshold & it does not drain externally
		if(summary[i,7]<merge_th & summary[i,6]>0){
			delete=c(delete,i)
			bas1=summary[i,1]
			bas2=summary[i,6]

			#replace numbers in the subbasin matrix
			ilist=which(subbasin==bas1)
			subbasin[ilist]=bas2

			#replace numbers in the subbasin area  matrix
			ilistA=which(subbasinA==bas1)
			subbasinA[ilistA]=bas2

			#adjust the summary matrix for the downstream basin
			summary[which(summary[,1]==bas2),7]=summary[which(summary[,1]==bas2),7] + summary[i,7]

			#Change the downstream basin number for any upstream basins to downstream basin
			uplist=which(summary[,6]==bas1)
			summary[uplist,6]=bas2

			#print(paste(bas1 , "- downstream:", bas2, "Upstream:" ))
			#print(uplist)

		} #end if
	} # end for

	if(is.null(delete)==F){
		summary=summary[-delete,]
	}

} # if merge >0


#test=rivers
#test[test==0]=NA
#maskcol=colorRampPalette(c('black', 'black'))
#image.plot(subbasinA, zlim=c(35,50))
#image.plot((test*2), add=T, col=maskcol(2), legend=F)

output_list=list("segments"=subbasin, "subbasins"=subbasinA, "RiverMask"=rivers, "summary"=summary)

return(output_list)

}
