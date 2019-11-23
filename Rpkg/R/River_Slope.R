#' Apply minimum slope and secondary scaling to river cells
#' 
#' A funtion that will apply a minimum slope threshold for the primary flow direciton along river cells. This function can also limit secondary slopes out of river cells to zero. 
#'
#' @param direction Nx by ny matrix of flow directions following the following convention (1=down, 2=left, 3=up, 4=right)
#' @param slopex Nx by ny matrix of slopes in the x direction (Note: these should be face centered slopes as woudl be calcualted with SlopeCalcStan)
#' @param slopey Nx by ny matrix of slopes in the y direction (Note: these should be face centered slopes as woudl be calcualted with SlopeCalcStan)
#' @param minslope Threshold for slope adjustment. Any primary direciton slope for a river cell will be adjusted such that abs(slope)>=minslope
#' @param RiverMask Nx by Ny matrix indicating the mask of river cells that the min slope will be applied to.  
#' @param Remove.Sec Flag, if set to True any secondary outlfows on river cells will be set to zero. Defaults to False. 
#' @return This will generate the following outputs:
#'    - slopex - The adjusted nx by ny matrix of slopex values
#'    - slopey - The adjusted nx by ny matrix of slopex values
#'    - adj_mask - Nx by ny matrix indicating cells that were ajusted (1=adjusted, 0=not adjsuted). This will be all the cells where inslope/outslop<adj_th
#'    - SlopeOutlet - Nx by ny matrix of the outlet slope for every grid cell 
#'    - SlopeOutletNew - Nx by ny matrix of the outlet slope for every grid cell after processing

#' @export
RivSlope=function(direction, slopex, slopey, minslope, RiverMask, Remove.Sec=FALSE){

nx=dim(direction)[1]
ny=dim(direction)[2]

#Colums: (1)deltax, (2)deltay, (3) direction number assuming you are looking downstream
kd=matrix(0, nrow=4, ncol=3) 
kd[,1]=c(0,-1,0,0)
kd[,2]=c(-1,0,0,0)
kd[,3]=c(1, 2, 3, 4)  
 

#setup outputs
outdirslope=outdirslopeNew=adj_mask=matrix(0, nrow=nx, ncol=ny)
slopexNew=slopex
slopeyNew=slopey

## Loop over the domain ajusting slopes along river cells as needed
#only looping over internal cells
for(j in 2:(ny-1)){
  for(i in 2:(nx-1)){
    if(RiverMask[i,j]==1 ){
      #jj=10
      #i=test[jj,1]
      #j=test[jj,2]
      
      #list all of the outflow slopes _ order: down, left, up, right
      sec_out=c( max(slopey[i,(j-1)],0), max(slopex[(i-1),j],0), -min(slopey[i,j],0), -min(slopex[i,j],0))
      
      #print(paste(i,j,"Outs:", sec_out[1], sec_out[2], sec_out[3], sec_out[4], "DIR", direction[i,j]))
      #print(max(sec_out[-direction[i,j]]))
      
      #adjust the outlet slope
      if(is.na(direction[i,j])==F){
        sec_out[direction[i,j]]=0 #zero out the primary flow direction outflow so this is jsut a list of secondary outflows
        
        #Set the primary direciton slope to be >=minslope
        if(direction[i,j]==1 & abs(slopey[i,(j-1)])<minslope){
          slopeyNew[i,(j-1)]=sign(slopey[i,(j-1)])*minslope
          outdirslope[i,j]=slopey[i,(j-1)]
          outdirslopeNew[i,j]=slopeyNew[i,(j-1)]
          adj_mask[i,j]=1
        } else if(direction[i,j]==2 & abs(slopex[(i-1),j])<minslope){
              slopexNew[(i-1),j]=sign(slopex[(i-1),j])*minslope
              outdirslope[i,j]=slopex[(i-1),j]
              outdirslopeNew[i,j]=slopexNew[(i-1),j]
              adj_mask[i,j]=1
        } else if(direction[i,j]==3 & abs(slopey[i,j])<minslope){
              slopeyNew[i,j]=sign(slopey[i,j])*minslope
              outdirslopeNew[i,j]=slopeyNew[i,j]
              adj_mask[i,j]=1
        }else if(direction[i,j]==4 & abs(slopex[i,j])<minslope){
              slopexNew[i,j]=sign(slopex[i,j])*minslope
              outdirslopeNew[i,j]=slopexNew[i,j]
              adj_mask[i,j]=1
        }
        
        #if Remove.sec is TRUE set any secondary outflow  slopes to 0
        if(Remove.Sec==TRUE){
          #print(paste(i,j, "True", sec_out))
          if(max(sec_out)>0){
            if(sec_out[1]>0){
              slopeyNew[(i+kd[1,1]), (j+kd[1,2])]=0
              #adj_mask[(i+kd[1,1]), (j+kd[1,2])]=2
              adj_mask[i,j]=adj_mask[i,j]+1
              }
            if(sec_out[2]>0){
              slopexNew[(i+kd[2,1]), (j+kd[2,2])]=0
              #adj_mask[(i+kd[2,1]), (j+kd[2,2])]=2
              adj_mask[i,j]=adj_mask[i,j]+1
              }
            if(sec_out[3]>0){
              slopeyNew[(i+kd[3,1]), (j+kd[3,2])]=0
              #adj_mask[(i+kd[3,1]), (j+kd[3,2])]=2
              adj_mask[i,j]=adj_mask[i,j]+1
              }
            if(sec_out[4]>0){
              slopexNew[(i+kd[4,1]), (j+kd[4,2])]=0
              #adj_mask[(i+kd[4,1]), (j+kd[4,2])]=2
              adj_mask[i,j]=adj_mask[i,j]+1
              }
          }
        } #end if Remove.Sec
        
      } #end if directions
    } # end if on mask 
  } #end for i
} #end for j


output_list=list("slopex"=slopexNew, "slopey"=slopeyNew, 
                 "adj_mask"=adj_mask, "SlopeOutlet"=outdirslope, "SlopeOutletNew"=outdirslopeNew)

return(output_list)
}