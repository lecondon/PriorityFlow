GetBorder=function(mask_mat){
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

	#GetBorder - Function that reads in a mask
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
