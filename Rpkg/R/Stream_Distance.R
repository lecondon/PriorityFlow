#'Find the distance to th nearest stream point following drainag directions
#'
#'This Function uses a stream mask and a flow direction file to determine the overland 
#'flow distance from any point in the domain to its nearest stream neigbor following
#'the definined primary flow directions. 

#' @inheritParams CalcSubbasins
#' @inheritParams RiverSmooth
#' @param streammask Mask with a value of 1 for every stream cell and 0 for all non stream cells. Refer to the
#' \code{\link{CalcSubbasins}} function to derive this mask. 
#' @param domainmask Option mask of the domain area with a value of 1 for cells inside the domain and 0 for cells
#' outside the domain. If no maks is provided it will default to using the entire rectangular domain.
#' @export
#' 
StreamDist=function(direction, streammask, domainmask, d4=c(1,2,3,4)){
  nx=dim(direction)[1]
  ny=dim(direction)[2]
  
  ku=matrix(0, nrow=4, ncol=3) 
  ku[,1]=c(0,1,0,-1)
  ku[,2]=c(1,0,-1,0)
  ku[,3]=c(1, 2, 3, 4)  
  
  if(missing(domainmask)){domainmask=matrix(1, nrow=nx, ncol=ny)} #default to processing everything
  
  #renumber the directions to 1=down, 2=left, 3=up, 4=right if a different numbering scheme was used
  dir2=direction
  if(d4[1]!=1){dir2[which(direction==d4[1])]=1}
  if(d4[2]!=2){dir2[which(direction==d4[2])]=2}
  if(d4[3]!=3){dir2[which(direction==d4[3])]=3}
  if(d4[4]!=4){dir2[which(direction==d4[4])]=4}
  
  #start a queue with every cell on the streammask
  #NOTE- This is going to ignore cells that drain to the outside of the domain
  #could add a border option later potentially
  queue=which(streammask==1, arr.ind=T)
  #distance=distance in cells along the flow directions to any stream cell
  #streamx and streamy are the x,y indices of each cells closests stream cell respectively
  distance=matrix(NA, nrow=nx, ncol=ny)
  
  #initialize distance to 0 along the stream mask
  #initialize streamx and stream y to the stream cell indices  along the stream mask
  distance[which(streammask==1)]=0
  streamy=matrix(rep(1:ny,nx), ncol=ny, byrow=T)
  streamx=matrix(rep(1:nx,ny), ncol=ny, byrow=F)
  streamx[which(streammask==0)]=NA
  streamy[which(streammask==0)]=NA
  active=TRUE

  
  while(active==T){
    indx=queue[1,1]
    indy=queue[1,2]
  
   queuetemp=NULL
   #Loop over four directions check for non-stream neighbors pointing to this cell
    for(d in 1:4){
      tempx=indx+ku[d,1]
      tempy=indy+ku[d,2]
      #print(c(tempx,tempy, dir2[tempx,tempy]))
      #if its point to the cell, is within the mask of cells to be processed, and has and epsilon < the threshold
      if(tempx*tempy>0 & tempx<nx & tempy<ny){
        if((d-dir2[tempx,tempy])==0 & streammask[tempx,tempy]==0 & domainmask[tempx,tempy]==1){
          #print(paste("CHECKING:", tempx,tempy, dir2[tempx,tempy],round(dem[tempx,tempy],1), round(dem[indx,indy],1)))
          distance[tempx,tempy]=distance[indx,indy]+1
          streamx[tempx,tempy]=streamx[indx,indy]
          streamy[tempx,tempy]=streamy[indx,indy]
          queuetemp=rbind(c(tempx,tempy),queuetemp)
        } #end drainge direction check
      } #end domain check
  }#end for d in 1:4
  
    #if cells were adjusted then add to the top of the queue replacing the cell that was just done
    if(length(queuetemp>0)){
      queue=rbind(queuetemp,queue[-1,])
    }else{
      #if no cells were adjsuted remove this cell from the queue and if its that last one you are done
      if(nrow(queue)>1){
        queue=queue[-1,]
        if(length(queue)==2){
          queue=matrix(queue, ncol=2, byrow=T)
        } #fixing bug to keep q formated as a mstrix if it drops down to one row
      }else{active=F}
    }
  } #end while active

  output_list=list("stream.dist"=distance, "stream.xind"=streamx, "stream.yind"=streamy)
  
  return(output_list)
}