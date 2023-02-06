#' Read PFB files
#' Function to Read ParFlow Binary Files into R
#' 
#' @param filename Name of file to read
#' @param verbose Optional to turn on print statement
#' @export 
readpfb = function(filename, verbose = FALSE){

  # Code to read parflow binary in R
  # RMM 10-19-13
  # JSAT 02-06-2023

  to.read = file(filename, "rb")
  on.exit(close(to.read))
  if(verbose){print(to.read)}
	
  # Read in X0, Y0, and Z0 of domain
  xyz = readBin(to.read, what = "double", n = 3, endian = "big")
  if (verbose){names(xyz) = c("X0", "Y0", "Z0"); print(xyz)}

  # Read in global NX, NY, and NZ of domain
  nxnynz = readBin(to.read, what = "integer", n = 3, endian = "big")
  if (verbose){names(nxnynz) = c("NX", "NY", "NZ"); print(nxnynz)}

  # Set up a blank array with the size of the domain
  data = array(0, dim = c(nxnynz[1], nxnynz[2], nxnynz[3]))

  # Read in DX, DY, and DZ
  dxdydz = readBin(to.read, what = "double", n = 3, endian = "big")
  if (verbose){names(dxdydz) = c("DX", "DY", "DZ"); print(dxdydz)}
	
  # Read in number of subgrids
  nsgs = readBin(to.read, what = "integer", n = 1, endian = "big")
  if (verbose){names(nsgs) = c("SUBGRIDS"); print(nsgs)}
  
  for (sg in 1:nsgs) {
    
    # Read in subgrid parameters: local starting points, local NX, NY, and NZ, subgrid refinement
    sgpar = readBin(to.read, what = "integer", n = 9, endian = "big")
    if (verbose){names(sgpar) = c("IX", "IY", "IZ", "NNX", "NNY", "NNZ", "RX", "RY", "RZ"); print(sgpar)}
    nr = sgpar[4] * sgpar[5] * sgpar[6]
    
    # Read in subgrid data
    x_inds = (sgpar[1] + 1):(sgpar[1] + sgpar[4])
    y_inds = (sgpar[2] + 1):(sgpar[2] + sgpar[5])
    z_inds = (sgpar[3] + 1):(sgpar[3] + sgpar[6])
    inds = as.matrix(expand.grid(x_inds, y_inds, z_inds))
    data[inds] = readBin(to.read, what = "double", n = nr, endian="big")
    
  }
  return(data)

}
