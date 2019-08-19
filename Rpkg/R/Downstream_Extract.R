#'Walk downstream from a point and extract values from a matrix
#'
#'This function grabs out values from a matrix by walking downstream from a point

#' @inheritParams CalcSubbasins
#' @param  input The matrix of values that you would like to extract the streampath from
#' @param  startpoint The x,y index of the grid cell you would like to start from 
#' @export
#' 
PathExtract=function(input, direction, mask, startpoint, d4=c(1,2,3,4)){
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
  nx=dim(direction)[1]
  ny=dim(direction)[2]
  
  if(missing(mask)){mask=matrix(1, nrow=nx, ncol=ny)} #default to processing everything

  path.mask=matrix(0, nrow=nx, ncol=ny) #matrix mapping the path
  path=NULL #List of cells on the path in order
  
  #D4 neighbors
  #Rows: down, left top right
  #Colums: (1)deltax, (2)deltay, direction number if you are waking (3)downstream, and  (4)upstream
  kd=matrix(0, nrow=4, ncol=4) 
  kd[,1]=c(0,-1,0,1)
  kd[,2]=c(-1,0,1,0)
  kd[,3]=c(1, 2, 3, 4)  
  kd[,3]=c(3, 3, 3, 2)  
  
  #renumber the directions to 1=down, 2=left, 3=up, 4=right if a different numbering scheme was used
  dir2=direction
  if(d4[1]!=1){dir2[which(direction==d4[1])]=1}
  if(d4[2]!=2){dir2[which(direction==d4[2])]=2}
  if(d4[3]!=3){dir2[which(direction==d4[3])]=3}
  if(d4[4]!=4){dir2[which(direction==d4[4])]=4}
  
  #intializing things
  indx=startpoint[1]
  indy=startpoint[2]
  step=1
  active=T
  output=NULL
  
  #walking downstream
  while(active==T){
    #print(paste(indx,indy,step))
    output=c(output,input[indx,indy])
    path.mask[indx,indy]=step
    path=rbind(path, c(indx,indy))
    
    #lookdownstream
    dirtemp=dir2[indx,indy]
    downindx=indx+kd[dirtemp,1]
    downindy=indy+kd[dirtemp,2]
    
    #If you have made it out of the domain then stop
    if(downindx<1 | downindx>nx | downindy<1 | downindy>ny){
      active=F
    }else{
      if(mask[downindx,downindy]==0){
        active=F
      }
    }
    
    #update the new indices
    indx=downindx
    indy=downindy
    step=step+1
  }
  
  output_list=list("data"=output, "path.mask"=path.mask, "path.list"=path)
  
}
  