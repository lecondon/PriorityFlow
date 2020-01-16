#' Initilze queue for topographic processing
#' 
#' Sets up a queue and initilizes marked and step matrices for DEM processing
#'
#' @inheritParams D4TraverseB
#' @param initmask Mask of the same dimensions of dem denoting a subset of cells to be considered for the queue (e.g. if you want to setup a run starting with only river cells). 
#' Note: if no init mask is included every border cell will be added to the queue
#' @param domainmask Mask of the domain extent to be considered. If no domain mask is provided boundaries will be calcualted from the rectangular extent
#' @param border Alternatively you can input your own border rathern than having it be calculated from the domain mask. For example if you want to have the river network and the borders combined you can input this as a border.
#' @return This function returns three outputs: "basins"=basin, "direction"=direction
#' @return 1.marked - A matrix indicating the outlet cells that were identified (1=outlet, 0=not outlet)
#' @return 2.queue - A list of the outlet cells. This has three columns: x, y, elevation 
#' This mask should match the river reach mask you provide the function if all cells were appropriately processed. 
#' @return 3.initmask- A matrix indicating the cells that were input as potential output points. 
#' @return 4.basin- A matrix indicating the basin number for each outlet point (each outlet is assigned a unique basin number)
#' @return 5.direciton- A matrix indicating the flow direction for each outlet point. The numbering scheme will follow the d4 numbering scheme provided as input.
#' @export
RunPF=function(dem, initmask, domainmask, border, d4=c(1,2,3,4)){
  ##Traverse the stream network
  t0=proc.time()
  trav1 = StreamTraverse(dem=demT, mask=channelmT, queue=init$queue, marked=init$marked, basins=init$basins, printstep=F, epsilon=ep)
  t1=proc.time()
  print(paste("First Pass:", round(t0[3]-t1[3],1), "sec"))
  print(paste("First Pass:", round(sum(trav1$marked)/sum(channelmT)*100,1), " % cells processed"))
  
  ##Look for orphaned branches and continue traversing until they are all connected
  norphan=1
  lap=1
  while(norphan>0){
    t1=proc.time()
    #look for orphan branches
    RivBorder=borderT+LborderT +SborderT+trav1$marked #mask of marked rivers+boundaries+lakes+sinks
    RivBorder[RivBorder>1]=1
    orphan=FindOrphan(trav1$dem, mask=channelmT, marked=RivBorder)
    norphan=orphan$norphan
    print(paste("lap", lap, norphan, "orphans found"))
    
    #go around again if orphans are found
    if(norphan>0){
      trav2 = StreamTraverse(trav1$dem, mask=channelmT, queue=orphan$queue, marked=trav1$marked, basins=trav1$basins, step=trav1$step, direction=trav1$direction, printstep=F, epsilon=ep)
      trav1=trav2
      lap=lap+1
      t2=proc.time()
      print(paste("Lap", lap , round(t2[3]-t1[3],1), "sec"))
    } else {
      print("Done!  No orphan branches found")
    }
    
  }
  print(paste("Final pass:", round(sum(trav1$marked * channelmT)/sum(channelmT)*100,1), " % cells processed"))
  t3=proc.time()
  print(paste("Total Time:" , round(t3[3]-t0[3],1), "sec"))
  
  ##Initialize the queue with every cell on the processed river and the boundary. 
  ###River border equals to the traversed river plus domain border plus lake and sink border.
  RivBorder=borderT+trav1$marked+ LborderT + SborderT 
  
  ###Initilize the updated river border.
  init=InitQueue(trav1$dem,  border=RivBorder)
  
  ## Process all the cells outside the channel network
  t4=proc.time()
  travHS=D4TraverseB(trav1$dem, queue=init$queue, marked=init$marked, mask=LakemaskT, direction=trav1$direction, basins=trav1$basins, step=trav1$step, epsilon=0.1, printstep=F, nchunk=1000) #adding lakes in the mask
  t5=proc.time()
  print(paste("Total Time:" , round(t5[3]-t4[3],1), "sec"))
}