write.raster=function(data, fout, xllcorner=0.0, yllcorner=0.0, dx=1000, naval=-999){
# This function writes out a matrix with the six header rows needed for GIS rasters
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
