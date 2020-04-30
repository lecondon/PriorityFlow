#' Defining Stream order
#' 
#' Function to calculate Strahler stream orers using stream segments delineated using the CalcSubbasins Function
#' @param basinID a column list of basin IDs
#' @param dsID downstream basin IDs, the downstream ID of each basin, should be correpsonding to the basin ID
#' @param segments (optional) a channel mask with segments assigned with basin IDs 
#' if PriorityFlow is used to generate subbasins, all inputs can be obtained with function CalcSubbasins in PriporityFlow.
#' 
#' @return summary - A summary table with a row for every basin with three columns "Basin_ID" ""Downstream_ID" "StreamOrder_number"
#' @return channel_orders (optional)- only available if segments is an input,this will output a spatial map with the channel orders
#' @export
CalStreamOrder=function(basinID,dsID,segments){

  order_no=array(0,c(length(basinID),1))
  
  #find 1st order, streams without any basin draining into this basin
  dsall=unique(dsID)
  hdwater=basinID[!basinID %in% dsall]
  for (i in 1:length(hdwater)){
    blist=which(basinID==hdwater[i])
    order_no[blist]=1
  }

  
  for (i in 1:length(hdwater)){
    active=T
    btemp=hdwater[i]

    while(active==T){
      #find downstream basin
      blist=which(basinID==btemp)
      dstemp=dsID[blist]
      dlist=which(basinID==dstemp)
      
      if(dstemp==0){active=F} #stop when the basin drains outside the domain

      #find all basins draining to this downstream basin
      ulist=which(dsID==dstemp)

      if(length(ulist)!=1){#if more than one basin
        urest=ulist[!ulist %in% blist]#remove the basin in process
        ordertemp=order_no[urest]
        ordertemp[is.na(ordertemp)]=0
        
        if (prod(ordertemp)!=0){#check if there is any upstreams not been processed
          umax=max(ordertemp)
          
          if (umax==order_no[blist]){#compare the orders of all tributaries
            order_no[dlist]=umax+1
          }else{
            order_no[dlist]=max(umax,order_no[blist])
          }
          
          btemp=dstemp
        }else{
          active=F
        }
      }else{#if only one is draining, assign the same order 
        order_no[dlist]=order_no[blist]
        btemp=dstemp
      }
    }
    
  }

  summary=cbind(basinID,dsID,order_no)
  colnames(summary)=c("Basin_ID","Downstream_ID","Order_number")
  
  if (missing(segments)){
    outputlist=list("summary"=summary)
  }else{
    segments2=segments
    for (i in 1:length(basinID)){
      btemp=basinID[i]
      blist2=which(segments2==btemp)
      segments[blist2]=order_no[i]
    }
    outputlist=list("summary"=summary, "channel_orders"=segments)
  }
  
  return(outputlist)
}

