#'Walk upstream from a point ensuring DEM is increasing by a minum epsilon
#'
#'This walks upstream from a point following a flow direction file and checks that the elevation upstream
#'is greater than or equal to the elvation at the point + epsilon
#'Once the function reaches a point where upstream cells pass it stops 
#'NOTE: this is only processing the immediate neigborhood and does not recurse over the entire domain, 
#'For example if you did the entire overal prioirty flow processing with an espilon of 0 then ran this
#'functon starting at a river point with an epsilon of 0.1 this function will traverse upstream from the 
#'river bottom checking that every cell conneting to the river cell is higher by at least this ammout but 
#'once it reaches a point where every connected cell passes this test it will stop. Therefore there could
#'still be locations higher up on the hillslope with the original epsilon of zero still.  This  is on purpose
#'and this script is not intended to gobally ensure a given epsilon. 

#' @inheritParams CalcSubbasins
#' @inheritParams RiverSmooth
#' @param startpoint the x,y index of a grid cell to start walking upstream from
#' @export
#' 
FixDrainage=function(dem, direction, mask, epsilon, startpoint, d4=c(1,2,3,4)){
  #D4 neighbors
  #Rows: down, left top right
  #Colums: (1)deltax, (2)deltay, direction number if you are waking upstream
  ku=matrix(0, nrow=4, ncol=3) 
  ku[,1]=c(0,1,0,-1)
  ku[,2]=c(1,0,-1,0)
  ku[,3]=c(1, 2, 3, 4)  
  
  #renumber the directions to 1=down, 2=left, 3=up, 4=right if a different numbering scheme was used
  dir2=direction
  nx=dim(direction)[1]
  ny=dim(direction)[2]
  if(d4[1]!=1){dir2[which(direction==d4[1])]=1}
  if(d4[2]!=2){dir2[which(direction==d4[2])]=2}
  if(d4[3]!=3){dir2[which(direction==d4[3])]=3}
  if(d4[4]!=4){dir2[which(direction==d4[4])]=4}
  
  #initializing 
  marked=matrix(0, nrow=nx, ncol=ny)
  queue=cbind(startpoint[1],startpoint[2])
  active=TRUE
  dem2=dem
  
  while(active==T){
    indx=queue[1,1]
    indy=queue[1,2]
    #print(paste("HERE", indx, indy))
    queuetemp=NULL
    #Loop over four directions check for non-stream neighbors pointing to this cell
    for(d in 1:4){
      tempx=indx+ku[d,1]
      tempy=indy+ku[d,2]
      #print(c(tempx,tempy, dir2[tempx,tempy]))
      #if its point to the cell, is within the mask of cells to be processed, and has and epsilon < the threshold
      if(tempx*tempy>0 & tempx<nx & tempy<ny){
        if((d-dir2[tempx,tempy])==0 & mask[tempx,tempy]==1){
            #print(paste("CHECKING:", tempx,tempy, dir2[tempx,tempy],round(dem[tempx,tempy],1), round(dem[indx,indy],1)))
           if((dem2[tempx,tempy]-dem2[indx,indy])<epsilon){
              dem2[tempx,tempy]=dem2[indx,indy]+epsilon
              #print(paste("ADJUSTING:", tempx,tempy, "FROM", round(dem[tempx,tempy],1), "TO:", round(dem2[tempx,tempy],1)))
              marked[tempx,tempy]=1
              queuetemp=rbind(c(tempx,tempy),queuetemp)
           } #end epsilon check
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
  #print(paste("Done! Start pont: ", startpoint[1], startpoint[2], ",", sum(marked), " cells adjusted", epsilon, sep=""))
  
  output_list=list("dem.adj"=dem2, "processed"=marked)
  
  return(output_list)

}
  
