FixBorderDir=function(direction, dem, d4=c(1,2,3,4)){
	#Oringinally all border cells are initialized to point out of the domain. This script checks the slopes and areas to determine if they need to swap


#If its a rectangle just calcualte the slopes the internal edges of the boundary cell and use this to determine which way it should go
	ny=ncol(direction)
	nx=nrow(direction)
		
	#Check 
	#top 
	#Check the slop between ny and ny-1 if it is positive then the flow should point in 
	sy_top=dem[,ny]-dem[,(ny-1)]
	flip=which(sy_top>0)
	direction[flip,ny]=d4[1]
	
	#bottom
	sy_bot=dem[,2]-dem[,1]
	flip=which(sy_bot<0)
	direction[flip,1]=d4[3]
	
	#right
	sx_right=dem[nx,]-dem[(nx-1),]
	flip=which(sx_right>0)
	direction[nx,flip]=d4[2]
	
	#left
	sx_left=dem[2,]-dem[1,]
	flip=which(sx_left<0)
	direction[1,flip]=d4[4]
	
	return(direction)
}