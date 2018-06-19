GetBorder=function(mask_mat){
	#This script reads in a mask
	#and returns a mask of the border cells for the irrigular boundary
	#The input mask should have values of 0 for cells outside the domain 1 for cells inside the domain
	
	ny=nrow(mask_mat)
	nx=ncol(mask_mat)
	border=matrix(1, nrow=ny, ncol=nx)
	border[2:(ny-1), 2:(nx-1)]= mask_mat[1:(ny-2), 2:(nx-1)] + 																mask_mat[3:ny, 2:(nx-1)] +
								mask_mat[2:(ny-1), 1:(nx-2)] +
								mask_mat[2:(ny-1), 3:nx]
	border=border*mask_mat
	border[which(border<4 & border!=0)]=1 
	border[border==4]=0
	
	return(border)
}