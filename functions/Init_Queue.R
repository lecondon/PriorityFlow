InitQueue=function(dem, initmask, domainmask, border){
	
#This function sets up a queue and initilizes marked and step
# Required Input: dem matrix transformed so that idices are [x,y]
# Optional Input:
# initmask - Mask of the same dimesions of dem denoting a subset of cells to be considered for the queue (e.g. if you want to setup a run starting with only river cells) 
# Note: if no mask is included every border cell will be added to the queue
# domainmask - mask of the domain extent to be considered. If no domain mask is 
# provided boundaries will be calcualted from the rectangular extent
# border- Alternatively you can input your own border rathern than having it be calculated from the domain mask. For example if you want to have the river network and the borders combined you can input this as a border. 

#initilize queue and matrices
ny=ncol(dem)
nx=nrow(dem)
queue=NULL
marked=matrix(0, nrow=nx, ncol=ny)
step=matrix(0, nrow=nx, ncol=ny)

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
	border[2:(nx-1), 2:(ny-1)]= domainmask[1:(nx-2), 2:(ny-1)] + 													domainmask[3:nx, 2:(ny-1)] +
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
	
output_list=list("mask"=initmask,"queue" = queue, "marked"=marked, "basins"=basin)
return(output_list)
}