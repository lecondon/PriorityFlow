#' Calculate Overland Flow
#' 
#' Function to calculate overland flow from ParFlow pressure file outputs. This function will write three pfb outputs 
#' for every timestep (1) outflow - this is the volumetric outflow from each grid cell [l3/t], (2) q_east - the 
#' volumetric flow across the east face of each cell (Note: this will have dimensions nx+1, ny to account for the western edge
#' of the domain) (2) q_north - the volumetric flow across the north face of each cell (Note this will have dimensions nx, ny+1
#'  to account for the southern edge of the domain)
#' @param file.path path where pressure and slope files are located, this is also where flow files will be written to
#' @param run.name ParFlow run name this is used to read the pressure files assumign they are named with the
#'         standard ParFlow naming convention (i.e. runname.out.press.00000.pfb)
#' @param file.nums list of file numbers to be processed
#' @param slopex.file name of the x slopes file  (note this must be a pfb file)
#' @param slopey.file  name of the y slopes file  (note this must be a pfb file)
#' @param overland.type  type of overland flow used for ParFlow simulation. Choices are 
#'        "OverlandFlow" or "OverlandKinematic" this shoudl match what you used in your tcl script. 
#' @param mannings manning roughness coefficient, can be either a value or a matrix
#' @param epsilon Epsilon used only in the OverlandKinematic formulation, this is set with 
#' the Solver.OverlandKinematic.Epsilon key in ParFlow (default value 1E-5)
#' @param dx grid size in x direciton [l] (default value 1)
#' @param dy grid size in y direciton [l] (default value 1)
#' @param  mask mask with ones for cells to be processed and zeros for everything else - defaults to a mask of all 1's
#'       if a mask is a provide it shoudl be a [nx,ny] matrix such that such that [1,1] is the lower left corner of the 
#'       domain and [nx,ny] is the upper right
#'       
#'       
CalcFlow=function(file.path, run.name, file.nums, slopex.file,
                  slopey.file, overland.type, mannings, epsilon=1e-5,
                  dx=1, dy=1, mask){
  
  
#Get the range of file number to be read in   
nfile=length(file.nums)

#make a zeros matrix for pmax calcs later
zeros=matrix(0, nrow=nx, ncol=ny)

#read in the slope files
fin.x=sprintf("%s/%s", file.path, slopex.file)
slopex=readpfb(fin.x, verbose=F)
fin.y=sprintf("%s/%s", file.path, slopey.file)
slopey=readpfb(fin.y, verbose=F)
slopex=slopex[,,1]
slopey=slopey[,,1]
nx=dim(slopex)[1]
ny=dim(slopex)[2]

if(missing(mask)){mask=matrix(1, nrow=nx, ncol=ny)} #default to processing everything
if(abs(nrow(mask)-nx)+abs(ncol(mask)-ny)!=0){
  stop("ERROR: X Y dimensions of slope and mask do not match")
}

for(f in 1:nfile){
  fn=file.nums[f]
  #read in the pressure file and get grid dimensions
  press.file=sprintf("%s/%s.out.press.%05d.pfb", file.path, run.name, fn)
  press=readpfb(press.file, verbose=F)
  nz=dim(press)[3]
  
  #Test that pressure files and slope files are the same size
  nxtest=dim(press)[1]
  nytest=dim(press)[1]
  if(abs(nxtest-nx)+abs(nytest-ny)!=0){
    stop("ERROR: X Y dimensions of slope and pressure files do not match")
  }

  #Get positive pressures in the top layer for overlandflow
  ptop=press[,,nz]
  ptop[ptop<0]=0 

  ################################################
  ## Overland Flow
  ################################################
  if(overland.type=="OverlandFlow"){
    #####
    #Calculate fluxes across east and north faces
    #First the x direction
    qx= -sign(slopex)*abs(slopex)^0.5/mannings * ptop^(5/3) *dy #Units should be l3/t
  
    #Upwinding to get flux across the east face of cells - based in qx[i] if its positive and qx[i+1] if its negative
    qeast= pmax(qx[1:(nx-1),],zeros[1:(nx-1),]) - pmax(-qx[2:nx,],zeros[2:nx,])
  
    #adding the left boundary - pressures outside domain are 0 so flux across this boundary only occurs when qx[1] is negative
    qeast= rbind(-pmax(-qx[1,],0), qeast)
  
    #adding the right boundary - pressures outside domain are 0 so flux across this boundary only occurs when qx[nx] is positive
    qeast= rbind(qeast, pmax(qx[nx,],0))
  
    #####
    #Next the y direction
    qy= -sign(slopey)*abs(slopey)^0.5/mannings * ptop^(5/3) * dx #Units should be l3/t
  
    #Upwinding to get flux across the north face of cells - based in qy[j] if its positive and qy[j+1] if its negative
    qnorth= pmax(qy[,1:(ny-1)],zeros[,1:(ny-1)]) - pmax(-qy[, 2:ny],zeros[, 2:ny])
  
    #adding the bottom - pressures outside domain are 0 so flux across this boundary only occurs when qy[1] is negative
    qnorth= cbind(-pmax(-qy[,1],0), qnorth)
  
    #adding the right boundary - pressures outside domain are 0 so flux across this boundary only occurs when qx[nx] is positive
    qnorth= cbind(qnorth, pmax(qy[,ny],0))
  
    #t(qeast[,ny:1])
    #t(qnorth[,(ny+1):1])
    #t(outflow[,ny:1])

  ################################################
  ## Overland Kinematic
  ################################################
  } else if (overland.type=="OverlandKinematic"){
    ####
    #Repeat the slopes on the lower and left boundaries that are inside the domain but outside the mask
    #find indices of all cells that are off the mask but have a neigbor to their right that is on the mask
    fill.left=which((rbind(mask[2:nx,],rep(0,ny)) - mask[1:nx,]) ==1, arr.ind=T)
    #get the indices of their neigbors to the right 
    fill.left2=fill.left
    fill.left2[,1]=fill.left[,1]+1
    #pad the slopes to the left with their neigboring cells in the mask 
    slopex[fill.left]=slopex[fill.left2]
    
    #find indices of all cells that are off the mask but have a neigbor above them that is on the mask
    fill.down=which((cbind(mask[,2:ny],rep(0,nx)) - mask[,1:ny]) ==1, arr.ind=T)
    #get the indices of their neigbors above
    fill.down2=fill.down
    fill.down2[,2]=fill.down[,2]+1
    #pad the slopes to below  with their neigboring cells in the mask 
    slopey[fill.down]=slopex[fill.down2]
    
    
    ####
    #calculate the slope magnitude
    sfmag=pmax((slopex^2+slopey^2)^0.5, epsilon)
    
    
    ###
    #For OverlandKinematic slopes are face centered and calculated across the upper and right boundaries
    # (i.e. Z[i+1]-Z[i])
    # For cells on the lower and left boundaries its assumed that the slopes repeat 
    # (i.e. repeating the upper and right face boundary for the lower and left for these border cells)
    slopex.pad=rbind(slopex[1,], slopex)
    slopey.pad=cbind(slopey[,1], slopey)
    
    ####
    # upwind the pressure - Note this is for the north and east face of all cells
    # The slopes are calculated across these boundaries so the upper boundary is included in these 
    # calculations and the lower and righ boundary of the domain will be added later
    pupwindx=pmax(sign(slopex)*rbind(ptop[2:(nx),], rep(0,ny)),0) +
           pmax(-sign(slopex)*ptop[1:nx,], 0 )
    pupwindy=pmax(sign(slopey) * cbind(ptop[,2:ny], rep(0,nx)),0) +
           pmax(-sign(slopey) * ptop[, 1:ny],0)
  
    ###
    # Calculate fluxes across east and north faces
    # First the x direction
    qeast = -slopex/(sfmag^0.5*mannings) * pupwindx^(5/3) *dy #Units should be l3/t
    qnorth = -slopey/(sfmag^0.5*mannings) * pupwindy^(5/3) *dx #Units should be l3/t
  
    ###
    #Fix the lower x boundary 
    # Use the slopes of the first column with the pressures for cell i 
    qleft=-slopex[1,]/(sfmag[1,]^0.5*mannings)* (pmax(sign(slopex[1,])*ptop[1,],0))^(5/3) * dy
    qeast=rbind(qleft,qeast)
  
    ###
    # Fix the lower y boundary 
    # Use the slopes of the bottom row with the pressures for cell j
    qbottom=-slopey[,1]/(sfmag[,1]^0.5*mannings)* (pmax(sign(slopey[,1])*ptop[,1],0))^(5/3) * dx
    qnorth=cbind(qbottom,qnorth)
  
} else {
  
  stop("Invalid overland.type: You must select either 'OverlandFlow' or 'OverlandKinematic'")
}


  #Calculate total outflow
  #Outflow is a postitive qeast[i,j] or qnorth[i,j] or a negative qeast[i-1,j], qnorth[i,j-1]
  outflow=pmax(qeast[2:(nx+1),],zeros) + pmax(-qeast[1:nx,], zeros) + 
          pmax(qnorth[,2:(ny+1)],zeros) + pmax(-qnorth[, 1:ny], zeros)

  #Write the outputs
  outflow.file=sprintf("%s/%s.out.outflow.%05d.pfb", file.path, run.name, fn)
  writepfb(outflow, outflow.file, dx, dy, 1)

  #Write the outputs
  qeast.file=sprintf("%s/%s.out.q_east.%05d.pfb", file.path, run.name, fn)
  writepfb(qeast, qeast.file, dx, dy, 1)

  #Write the outputs
  qnorth.file=sprintf("%s/%s.out.q_east.%05d.pfb", file.path, run.name, fn)
  writepfb(qnorth, qnorth.file, dx, dy, 1)

} #end for f
#output_list=list("qeast"=qeast, "qnorth"=qnorth, "outflow"=outflow)
#return(output_list)

}