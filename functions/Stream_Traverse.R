StreamTraverse=function(dem, mask, queue, marked, step, direction,basins, d4=c(1,2,3,4), printstep=F, epsilon=0){

#Function to process stream networks walking upstream on d4 neigbors
#in a river mask. Where no D4 neigbhors exist it looks for d8 neigbors
#and created d4 bridges to these diagonal cells

#Mandatory Inputs:
# 1. dem: Elevation matrix
# 2. mask: Mask with zeros for non river cells and 1 for river cells
# 3. queue: a priority queue to start from three columns, x, y, elevation
# 4. marked: a matrix of which cells have been marked already

#Optional Inputs:
# 1. d4: directional numbering system: the numbers 
#	you want to assigne to down, left, top,right
#   defaults to 1,2,3,4
# 2. printstep: if true it will print out the step number and the size of the queue
# 3. epsilon: amount to add to filled areas to avoid creating flats, defaults to zero 
# 4. step: a matrix of the step number for cells that have been processed - defaults to all zeros
# 5. direction: a matrix of hte flow directions for cells that have been processed - defaults to all zeros
# 6. basins: a matrix of basin numbers that can be created by the initilizaiton script. If you input this every cell will be assigned the same basin as the cell that adds it



nx=dim(dem)[1]
ny=dim(dem)[2]
demnew=dem

#setup matrices for anything that wasn't input
if(missing(step)){step=matrix(0, nrow=nx, ncol=ny)} #start all steps at zero
if(missing(direction)){direction=matrix(NA,nrow=nx, ncol=ny)} #make a blank direction matrix
if(missing(basins)){basins=matrix(0,nrow=nx, ncol=ny)} #make all the basins=1

#D4 neighbors
kd=matrix(0, nrow=4, ncol=3) #ordered down, left top right 
kd[,1]=c(0,-1,0,1)
kd[,2]=c(-1,0,1,0)
kd[,3]=c(d4[3], d4[4], d4[1], d4[2])  #We are walking upstream so the direction needs to point opposite

#D8 neighbors
kd8=matrix(0, nrow=4, ncol=4)
kd8[,1]=c(-1,-1,1,1)
kd8[,2]=c(-1,1,1,-1)
#directions for the two d4 neighbor cells that are tested
#the first neigbhor cell is xC,yk (i.e. just moving up or down)
kd8[,3]= c(d4[3], d4[1], d4[1], d4[3]) 
#the second neigbor dell is xk, yC (i.e. just moving left or right)  
kd8[,4]= c(d4[4], d4[4], d4[2], d4[2]) 


nqueue=nrow(queue)
nstep=0

while(nqueue>0){

##############
#pick the lowest DEM cell on the queue
pick=which.min(queue[,3])
xC=queue[pick,1]
yC=queue[pick,2]
demC=queue[pick,3]

#print(paste('pick:', xC, yC))


##############
# Look for D4 neighbor cells that are on the mask and add to queue
count=0
for(k in 1:4){
	xk=xC+kd[k,1]
	yk=yC+kd[k,2]
	#print(c(xk, yk))
	#check that the neigbhor is inside the domain and on the mask
	if(yk>=1 & yk<=ny & xk>=1 & xk<=nx ){
			if (mask[xk, yk]==1 & marked[xk,yk]==0){
				#print(c(k, xk,yk))
				demtemp=max((demC+epsilon), dem[xk, yk])
				demnew[xk, yk]=demtemp
				queue=rbind(queue, c(xk, yk, demtemp))
				marked[xk, yk]=1
				step[xk, yk]=step[xC, yC]+1
				direction[xk,yk]=kd[k,3]
				basins[xk,yk]=basins[xC,yC]
				count=count+1
				nqueue=nqueue+1
			}
		}		
}
#print(paste(count, "available D4 neighbors found"))

##############
# If you don't find any D4 neigbors look for D8 neigbors
# and choose the least cost D4 option to reach that D8 cell
n4=matrix(NA,ncol=4, nrow=2)
if(count==0){
	#print("No D4 neighbors, checking for diagonal")
	for(k in 1:4){
		xk=xC+kd8[k,1]
		yk=yC+kd8[k,2]
		#print(c(xk, yk))
		count4=0
		n4=matrix(NA,ncol=4, nrow=2)
		#check that the neigbhor is inside the domain and on the mask
		if(yk>=1 & yk<=ny & xk>=1 & xk<=nx ){
			#check that it hasn't alredy been marked
			if(marked[xk,yk]<1 & mask[xk,yk]==1){
				#print("found 1!")
				#look for available d4 neigbhors to add instead
				if(marked[xC,yk]<1){
					n4[1,]=c(xC, yk, dem[xC, yk], kd8[k,3])
					count4=count4+1
				}
				if(marked[xk,yC]<1){
					n4[2,]=c(xk, yC, dem[xk, yC], kd8[k,4])
					count4=count4+1
				}
				
				#choose the neighbor which is the lowest without going under
				#and if not fill the highest
				if(count4>0){
					if(min(n4[,3], na.rm=T)>=demC){
						npick=which.min(n4[,3])
					} else{
						npick=which.max(n4[,3])	
					}
					demtemp=max((demC+epsilon), dem[n4[npick,1], n4[npick,2]])
					demnew[n4[npick,1], n4[npick,2]]=demtemp
					queue=rbind(queue, c(n4[npick,1], n4[npick,2], demtemp))
					marked[n4[npick,1], n4[npick,2]]=1
					direction[n4[npick,1], n4[npick,2]]=n4[npick,4]
					step[n4[npick,1], n4[npick,2]]=step[xC, yC]+1
					basins[n4[npick,1], n4[npick,2]]=basins[xC, yC]
					count=count+1
					#print(n4[npick,1:2])
				} else{
					#print(paste(xC, yC, 'D8 mask cell found with no viable d4 neighbors'))
				} #end if count4>0
				
			} #end if d8 is nt marked and on mask
		} #end if d8 is within the domain
	} #end for k in 1:4 (i.e. loping over d8 neighbors)
} # end if count = 0

##############
#Remove from the queue and move on
nqueue=length(queue)/3
if(nqueue>1){
	queue=queue[-pick,]
	nqueue=length(queue)/3
} else {
	nqueue=0
}

#if there is only one row it will treat it as a vector and things will break
if(nqueue==1){
	queue=matrix(queue, ncol=3,byrow=T)
}
 
nstep=nstep+1
if(printstep){print(paste("Step:", nstep, "NQueue:", nqueue))}

}



output_list=list("dem"=demnew, "mask"=mask, "marked"=marked, "step"= step, "direction"=direction, "basins"=basins)
return(output_list)

} # end function

