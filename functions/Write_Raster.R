write.raster=function(data, fout, xllcorner=0.0, yllcorner=0.0, dx=1000, naval=-999){

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

# Write.Raster - This function writes out a matrix with the six header rows needed for GIS rasters
# Required inputs:
# 	data - Matrix of values for the raster
#   fout - file name for output
# Optional inputs
# 	xllcorner - xlocation of lower left corner of domain (defaults to 0.0)
#   yllcorner - ylocation of lower left corner of domain (defaults to 0.0)
# 	naval - value assigned to NAs in the raster (defaults to -999)

header1=paste("ncols", ncol(data), sep="      ")
header2=paste("nrows", nrow(data), sep="      ")
header3=paste("xllcorner", xllcorner, sep="      ")
header4=paste("yllcorner", yllcorner, sep="      ")
header5=paste("cellsize", dx, sep="      ")
header6=paste("NODATA_value", naval, sep="      ")

header = rbind(header1, header2, header3, header4, header5, header6)

write.table(header, fout, row.names=F, col.names=F, append=F, quote=F)
write.table(data, fout, row.names=F, col.names=F, append=T)

}
