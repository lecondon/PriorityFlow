#' Indentify and fix stagnation points
#' 
#' A funtion that will apply a minimum slope threshold for the primary flow direciton along river cells. 
#'
#' @param direction Nx by ny matrix of flow directions following the following convention (1=down, 2=left, 3=up, 4=right)
#' @param slopex Nx by ny matrix of slopes in the x direction (Note: these should be face centered slopes as woudl be calcualted with SlopeCalcStan)
#' @param slopey Nx by ny matrix of slopes in the y direction (Note: these should be face centered slopes as woudl be calcualted with SlopeCalcStan)
#' @param minslope Threshold for slope adjustment. Any primary direciton slope for a river cell will be adjusted such that abs(slope)>=minslope
#' @param RiverMask Nx by Ny matrix indicating the mask of river cells that the min slope will be applied to.  
#' @return This will generate the following outputs:
#'    - slopex - The adjusted nx by ny matrix of slopex values
#'    - slopey - The adjusted nx by ny matrix of slopex values
#'    - adj_mask - Nx by ny matrix indicating cells that were ajusted (1=adjusted, 0=not adjsuted). This will be all the cells where inslope/outslop<adj_th
#'    - SlopeOutlet - Nx by ny matrix of the outlet slope for every grid cell 
#'    - SlopeOutletNew - Nx by ny matrix of the outlet slope for every grid cell after processing

#' @export
RivSlope=function(direction, slopex, slopey, minslope, RiverMask){

nx=dim(direction)[1]
ny=dim(direction)[2]

#setup outputs
outdirslope=outdirslopeNew=adj_mask=matrix(0, nrow=nx, ncol=ny)
slopexNew=slopex
slopeyNew=slopey

## Loop over the domain calculating the total input slope, output slope, and the slope at the outlet
#only looping over internal cells
for(j in 2:(ny-1)){
  for(i in 2:(nx-1)){
    if(RiverMask[i,j]==1 ){
     #adjust the outlet slope
      if(is.na(direction[i,j])==F){
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
      } #end if directions
    } # end if on mask 
  } #end for i
} #end for j


output_list=list("slopex"=slopexNew, "slopey"=slopeyNew, 
                 "adj_mask"=adj_mask, "SlopeOutlet"=outdirslope, "SlopeOutletNew"=outdirslopeNew)

return(output_list)
}