#' Indentify and fix stagnation points
#' 
#' A funtion that finds cells where the total output slope is much less than the input
#' slopes and adjusts the outlet slope to address this
#'
#' @param direction Nx by ny matrix of flow directions following the following convention (1=down, 2=left, 3=up, 4=right)
#' @param slopex Nx by ny matrix of slopes in the x direction (Note: these should be face centered slopes as woudl be calcualted with SlopeCalcStan)
#' @param slopey Nx by ny matrix of slopes in the y direction (Note: these should be face centered slopes as woudl be calcualted with SlopeCalcStan)
#' @param adj_th Threshold for slope adjustment. If the total slopes out divided by the total slopes into a given cell is less than adj_th the outlet will be scale by adj_ratio.
#' @param adj_ratio Scaler value for slope adjustment. New outlet slopes will be set to inital outlet slope time adj_ratio. If no ad_ratio is provided
#' then the Outlet slope will be set to the total inlets slopws times the adjustment theshold
#' @param mask Nx by Ny matrix indicating the active domain, 1=active, 0 in inactive. If not mask is provided the function will be applied over the entire input matrices. 
#' @return This will generate the following outputs:
#'    - outslope- Nx by ny matrix with the sum of the slopes pointed out of every cell (Note: this is the sum of slop magnitudes so it will all be positive)
#'    - inslope- Nx by ny matrix with the sum of the slopes pointed into every cell (Note: this is the sum of slop magnitudes so it will all be positive)
#'    - OutIn_ratio - Nx by Ny matrix with the ratio of outslope to inslope (Note: cells outside the mask or with 0 inslope are assigned a ratio of -1)
#'    - slopex - The adjusted nx by ny matrix of slopex values
#'    - slopey - The adjusted nx by ny matrix of slopex values
#'    - adj_mask - Nx by ny matrix indicating cells that were ajusted (1=adjusted, 0=not adjsuted). This will be all the cells where inslope/outslop<adj_th
#'    - SlopeOutlet - Nx by ny matrix of the outlet slope for every grid cell 
#'    - SlopeOutletNew - Nx by ny matrix of the outlet slope for every grid cell after processing

#' @export
FixFlat=function(direction, slopex, slopey, adj_th, adj_ratio=-1, mask){

nx=dim(direction)[1]
ny=dim(direction)[2]

if(missing(mask)){
  print("No domain mask provided using entire domain")
  mask=matrix(1, nrow=nx, ncol=ny)
}

#setup outputs
intot=outtot=outdirslop=ratio_mask=matrix(0, nrow=nx, ncol=ny)
outdirslope=outdirslopeNew=matrix(0, nrow=nx, ncol=ny)
ratio=matrix(-1, nrow=nx, ncol=ny) #setting background ratio value to -1
slopexNew=slopex
slopeyNew=slopey

## Loop over the domain calculating the total input slope, output slope, and the slope at the outlet
#only looping over internal cells
for(j in 2:(ny-1)){
  for(i in 2:(nx-1)){
    if(mask[i,j]==1){
      #total in slopes and out slopes
      intot[i,j]=max(slopex[i,j],0)+max(slopey[i,j],0)-min(slopex[(i-1),j],0)-min(slopey[i,(j-1)],0)
      outtot[i,j]=max(slopex[(i-1),j],0)+max(slopey[i,(j-1)],0)-min(slopex[i,j],0)-min(slopey[i,j],0)
    
      #ratio of out to in 
      ratio[i,j]=outtot[i,j]/intot[i,j]
      
      #adjust for cells where inflow is zero
      if(abs(ratio[i,j])==Inf){ratio[i,j]=(-1)}
      
      #if ratio is less than the adjustment threshold than flag it 
      if(ratio[i,j]<adj_th & ratio[i,j]>0){ratio_mask[i,j]=1}
      
     #adjust the outlet slope
      if(is.na(direction[i,j])==F){
        if(direction[i,j]==1){
          outdirslope[i,j]=slopey[i,(j-1)]
          #if the ratio is below the threshold then adjust outlet 
          if(ratio_mask[i,j]==1){
            if(adj_ratio>0){
              slopeyNew[i,(j-1)]=slopey[i,(j-1)]*adj_ratio
            } else{ 
              slopeyNew[i,(j-1)]=intot[i,j]*adj_th*sign(slopey[i,(j-1)])
            }
            outdirslopeNew[i,j]=slopeyNew[i,(j-1)]
          }
          
        } else if(direction[i,j]==2){
          outdirslope[i,j]=slopex[(i-1),j]
          #if the ratio is below the threshold then adjust outlet 
          if(ratio_mask[i,j]==1){
            if(adj_ratio>0){
              slopexNew[(i-1),j]=slopex[(i-1),j]*adj_ratio
            } else{ 
              slopexNew[(i-1),j]=intot[i,j]*adj_th*sign(slopex[(i-1),j])
            }
            outdirslopeNew[i,j]=slopexNew[(i-1),j]
          }
          
        } else if(direction[i,j]==3){
          outdirslope[i,j]=slopey[i,j]
          #if the ratio is below the threshold then adjust outlet 
          if(ratio_mask[i,j]==1){
            if(adj_ratio>0){
              slopeyNew[i,j]=slopey[i,j]*adj_ratio
            } else{ 
              slopeyNew[i,j]=intot[i,j]*adj_th*sign(slopey[i,j])
            }
            outdirslopeNew[i,j]=slopeyNew[i,j]
          }
          
        }else if(direction[i,j]==4){
          outdirslope[i,j]=slopex[i,j]
          #if the ratio is below the threshold then adjust outlet 
          if(ratio_mask[i,j]==1){
            if(adj_ratio>0){
              slopexNew[i,j]=slopex[i,j]*adj_ratio
            } else{ 
              slopexNew[i,j]=intot[i,j]*adj_th*sign(slopex[i,j])
            }
            outdirslopeNew[i,j]=slopexNew[i,j]
          }
        }
      } #end if directions
    } # end if on mask 
  } #end for i
} #end for j


output_list=list("outslope"=outtot, "inslope"=intot, "OutIn_ratio"=ratio, "slopex"=slopexNew, "slopey"=slopeyNew, 
                 "adj_mask"=ratio_mask, "SlopeOutlet"=outdirslope, "SlopeOutletNew"=outdirslopeNew)

return(output_list)
}