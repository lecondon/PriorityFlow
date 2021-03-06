#'Apply smoothing to a DEM along a pre-defined stream network
#'
#'This function will smooth a DEM along a stream network. It requires pre-defined stream segments 
#'and subbasins which can be obtained using the \code{\link{CalcSubbasins}} function.

#' @inheritParams CalcSubbasins
#' @inheritParams D4TraverseB
#' @param  river.summary A table summarizing the stream segments in a domain with the following. 
#'  it should have 1 row per stream segmentand the following 7 columns: (1) subbasin number, (2) x and (3) y index
#'  of the upstream end of the stream segment, (4) x and (5) y index of the downstream eand of the river segement,
#'  (6) the subbasin number for the downstream basin, -1 indicates a subbasin draining out of the domain
#'  (7) drainage area of the subbasin 
#' @param river.segments A nx by ny matrix indicating the subbasin number for with for all grid cells on the river network 
#'      ( (all cells not on the river network shoudl be 0)
#' @param bank.epsilon the minimum elevation difference between cells walking up the banks from the river network. 
#' @param river.epsilon the minimum elevation difference between cells along the river
#' @return This function returns three outputs:
#' @return 1.dem.adj - A matrix with the adjusted DEM values following the river smoothing opeartion
#' @return 2.processed - A matrix indicating the cells that were processed by this routine (1=processed, 0=not processed). 
#' This mask should match the river reach mask you provide the function if all cells were appropriately processed. 
#' @return 3.summary - A summary of the reach properties. This is a matrix with a row for every river reach and the following columns:
#' 1. River segment ID number, 2. X index of the top of the segment, 3. X index of the bottom of the segment, 
#' 4. Y index of the top of the segment, 5. Yindex of the bottom of the segment, 6. Length of the segment,
#' 7. Elevation at the top of the segment, 8. Elevation at the bottom of the segment, 9. Delta applied along the segment (i.e. (Top-Bottom)/Length)
#' @export
#' 
RiverSmooth=function(dem, direction, mask, river.summary, river.segments, bank.epsilon=0.01, river.epsilon=0.0,  d4=c(1,2,3,4), printflag=F){
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
  ##TEMP
  #dem=travHS$dem
  #river.summary=subbasin$summary
  #river.segments=subbasin$segments
  #subasins=subbasin$subbasins

  ##################
  #source("~/Dropbox/CONUS_Share/Topography_Testing/SmallTestDomain/Slope_processingNew/functions/Fix_Drainage.R", local=TRUE)
  #D4 neighbors
  #Rows: down, left top right
  #Colums: (1)deltax, (2)deltay, direction number if you are waking (3)downstream, and  (4)upstream
  kd=matrix(0, nrow=4, ncol=4) 
  kd[,1]=c(0,-1,0,1)
  kd[,2]=c(-1,0,1,0)
  kd[,3]=c(1, 2, 3, 4)  
  kd[,3]=c(3, 3, 3, 2)  
  
  if(missing(mask)){mask=matrix(1, nrow=nx, ncol=ny)} #default to processing everything

  #renumber the directions to 1=down, 2=left, 3=up, 4=right if a different numbering scheme was used
  dir2=direction
  nx=dim(direction)[1]
  ny=dim(direction)[2]
  if(d4[1]!=1){dir2[which(direction==d4[1])]=1}
  if(d4[2]!=2){dir2[which(direction==d4[2])]=2}
  if(d4[3]!=3){dir2[which(direction==d4[3])]=3}
  if(d4[4]!=4){dir2[which(direction==d4[4])]=4}
  
  #setup a river list 
  nriver=nrow(river.summary)
  marked.segments=rep(0,nriver) #marker for keeping track of which reaches are processed
  marked.matrix=matrix(0, ncol=ny, nrow=nx)
  
  #Setup a smoothing summary
  riversmooth.summary=matrix(0, nrow=nriver, ncol=9)
  colnames(riversmooth.summary)=c("SegmentID", "Start.X", "Start.Y", "End.X", "End.Y","Lenght", "Top.Elevation", "Bottom.Elevation", "delta")
  riversmooth.summary[,1:5]=river.summary[,1:5]
  
  #make a mask of the hillslope cells
  hillmask=mask
  hillmask[which(river.segments>0)]=0
  
  #First make a list of all the termial river reaches
  #queue=which(river.summary[,6]==(-1))
  queue=which(river.summary[,6]<=(0))
  if(length(queue)>0){active=TRUE}else{print("No terminal river segments provided, not adjusting DEM")}
  
  #start a new dem
  dem2=dem
  
  #get the length of every river segment
  river.length=rep(0,max(river.summary[,1])) #this has to be as long as the maximum river index not the number of rivers becasue there can be river indeces that don't exist in the river summary after merges occur
  for(i in 1:nx){
    for(j in 1:ny){
      rtemp=river.segments[i,j]
      river.length[rtemp]=river.length[rtemp]+1
    }
  }
  
  #Loop over the river segments working upstream
  #starting from every terminal segment (i.e. segments with a downstream river number of -1)
  while(active==T){
    indr=queue[1]
    r=river.summary[indr,1] #river segment number
    rdown=river.summary[indr,6]
    length=river.length[r]
    riversmooth.summary[indr,6]=river.length[r]
    
    #find the top and bottom elevations of the current river segment
    top=dem2[river.summary[indr,2], river.summary[indr,3]]
    #if its a terminal reach then the bottom elevation will be the bottom of the reach
    #if(rdown==(-1)){
    if(rdown<=0){
      bottom=dem2[river.summary[indr,4], river.summary[indr,5]]
      length=length-1
    } else{
    #if not then use the elevation downstream of the bottom point of the reach
      bdir=dir2[river.summary[indr,4], river.summary[indr,5]]
      bottom=dem2[(river.summary[indr,4]+kd[bdir,1]), (river.summary[indr,5]+kd[bdir,2])]
    }

    topmin=bottom+river.epsilon*length
    if(top<topmin){
      #calculate the delta from the original dem
      top0=dem[river.summary[indr,2], river.summary[indr,3]]
      if(rdown>0){
       bdir=dir2[river.summary[indr,4], river.summary[indr,5]]
        bottom0=dem[(river.summary[indr,4]+kd[bdir,1]), (river.summary[indr,5]+kd[bdir,2])]
      }else{
        bottom0=dem[river.summary[indr,4], river.summary[indr,5]]
      }
      
      #use this delta from the original dem to adjust the top elevation
      delta=max((top0-bottom0)/(length),river.epsilon)
      top=bottom+delta*length
      dem2[river.summary[indr,2], river.summary[indr,3]]=top
      if(printflag==T){
        print(paste("River top elevation<river bottom elevation for segment", r))
        print(paste("Original top", round(top0,2), "and original bottom", round(bottom0,2)))
        print(paste("Adjusting the top elevation from", round(top0,2), "top", round(top,2)))
      }
    }

    
    if(printflag==T){
      print(paste("River segment:",r))
      print(paste("Start:",river.summary[indr,2], river.summary[indr,3], round(top,1)))
      print(paste("End:",river.summary[indr,4], river.summary[indr,5], round(bottom,1)))
    }

    #walk from top to bottom smoothing out the river cells
    indx=river.summary[indr,2]
    indy=river.summary[indr,3]
    marked.matrix[indx,indy]=marked.matrix[indx,indy]+1

    if(length>1){delta=(top-bottom)/(length)}else{delta=0}
    if(delta<0){
      print(paste("Warning: Calculated delta < 0, setting delta to 0 for segment", r))
      delta=0
    }
    temp=top
    riversmooth.summary[indr,7]=top
    riversmooth.summary[indr,8]=bottom
    riversmooth.summary[indr,9]=delta
    
    if(length>1){
      for(i in 2:length){
        #print(i)
        temp=temp-delta
        #find the downstream point and adjust its elevation
        dirtemp=dir2[indx,indy]
        downindx=indx+kd[dirtemp,1]
        downindy=indy+kd[dirtemp,2]
        if(river.segments[downindx,downindy]==r){
          dem2[downindx,downindy]=temp
          marked.matrix[downindx,downindy]=marked.matrix[downindx,downindy]+1
          #loop up the hillslope from the point and make sure everything drains
          drainfix=FixDrainage(dem=dem2, direction=dir2, mask=hillmask, bank.epsilon=bank.epsilon, startpoint=c(downindx,downindy))
          dem2=drainfix$dem.adj
        }else{
          print(paste("Warnig: Check Segmennt for  branches", r))
        }
       
        #print(paste(downindx,  downindy, temp))
        #move to the downstream point
        indx=downindx
        indy=downindy
        
        

        #print(paste(indx, indy, round(temp)))
      } #end for i in length river segment
    } # end if length >1
    
    #mark this segement as done
    marked.segments[indr]= marked.segments[indr]+1
      
    #Find all of the river segments that drain to this segment
    uplist=which(river.summary[,6]==r)
    #if there are upstream segments then add to the top of the Q overwriting the segment that
    #was just done if not just remove the current segment from the queue
    if(length(uplist>0)){queue=c(uplist,queue[-1])} else{queue=queue[-1]}
    if(length(queue)==0){active=F}
  } #end while active
      
output_list=list("dem.adj"=dem2, "processed"=marked.matrix, "summary"=riversmooth.summary)
  
return(output_list)
} # End function
  
    
    
    
   

