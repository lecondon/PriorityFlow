D4TraverseB=function(dem, queue, marked, mask, step, direction, basins, d4=c(1,2,3,4), printstep=F, nchunk=100, epsilon=0){

#Function to process stream networks walking upstream on d4 neigbors
#in a river mask. Where no D4 neigbhors exist it looks for d8 neigbors
#and created d4 bridges to these diagonal cells

#Mandatory Inputs:
# 1. dem: Elevation matrix
# 2. queue: a priority queue to start from three columns, x, y, elevation
# 3. marked: a matrix of which cells have been marked already

#Optional Inputs:
# 1. d4: directional numbering system: the numbers 
#	you want to assigne to down, left, top,right
#   defaults to 1,2,3,4
# 2. printstep: if true it will print out the step number and the size of the queue
# 3. epsilon: amount to add to filled areas to avoid creating flats, defaults to zero 
# 4. mask: Mask with ones for cells to be processed and zeros for everything else - defaults to a mask of all 1's
# 5. step: a matrix of the step number for cells that have been processed - defaults to all zeros
# 6. direction: a matrix of the flow directions for cells that have been processed - defaults to all zeros
# 7. basins: a matrix of basin numbers that can be created by the initilizaiton script. If you input this every cell will be assigned the same basin as the cell that adds it

t0=proc.time()
nx=dim(dem)[1]
ny=dim(dem)[2]
demnew=dem

#setup matrices for anything that wasn't input
if(missing(mask)){mask=matrix(1, nrow=nx, ncol=ny)} #default to processing everything
if(missing(step)){step=matrix(0, nrow=nx, ncol=ny)} #start all steps at zero
if(missing(direction)){direction=matrix(NA,nrow=nx, ncol=ny)} #make a blank direction matrix
if(missing(basins)){basins=matrix(0,nrow=nx, ncol=ny)} #make all the basins=1

#D4 neighbors
kd=matrix(0, nrow=4, ncol=3) #ordered down, left top right 
kd[,1]=c(0,-1,0,1)
kd[,2]=c(-1,0,1,0)
kd[,3]=c(d4[3], d4[4], d4[1], d4[2])  #We are walking upstream so the direction needs to point opposite

split=0
q1max=0
nqueue=nrow(queue)
nstep=0
queuetemp=NULL

#split the queue in 2 using the top nchuck values for the first 
#queue and the rest for the second
if(nqueue>nchunk){
	qsort=queue[order(queue[,3]),]
	queue1=qsort[1:nchunk,]
	queue2=qsort[-(1:nchunk),]
	th=queue2[1,3]
	nqueue2=length(queue2)/3
	nqueue=nrow(queue1)
	print(paste('inital queue:', nrow(queue), "splitting. Q1=", nqueue, "Q2=", nqueue2))
	#adding the matrix step so it doesn't become a vector if its only one
	if(nqueue2==1){
		queue2=matrix(queue2, ncol=3,byrow=T)
	}
} else{
	print(paste('inital queue:', nqueue, "Not splitting"))
	queue1=queue
	nqueue2=0
}
t0=proc.time()

while(nqueue>0){

#for(jj in 1:2150){

#t0=proc.time()
##############
#pick the lowest DEM cell on the queue
#t1=proc.time()
pick=which.min(queue1[,3])
#t2=proc.time()
xC=queue1[pick,1]
yC=queue1[pick,2]
demC=queue1[pick,3]
#t1=proc.time()
#t1-t0

#print(paste('pick:', xC, yC))


##############
# Look for D4 neighbor cells that are on the mask and add to queue

count=0
for(k in 1:4){
	xk=xC+kd[k,1]
	yk=yC+kd[k,2]
	#print(c(xk, yk))
	#check that the neigbhor is inside the domain and on the mask
	#t3=proc.time()
	if(yk>=1 & yk<=ny & xk>=1 & xk<=nx ){
			if (mask[xk, yk]==1 & marked[xk,yk]==0){
				#print(c(k, xk,yk))
				demtemp=max((demC+epsilon), dem[xk, yk])
				demnew[xk, yk]=demtemp
				#if its less than the threshold add it to Q1
				# if not add it to Q2
				if(demtemp<th){
					queue1=rbind(queue1, c(xk, yk, demtemp))
					#print(paste("Q1 add:", round(demtemp,1), round(th,1)))
				} else{
					#queue2=rbind(queue2, c(xk, yk, demtemp))
					queuetemp=rbind(queuetemp, c(xk, yk, demtemp))
					#print(paste("Q2 add:", round(demtemp,1), round(th,1)))

				}
				
				marked[xk, yk]=1
				step[xk, yk]=step[xC, yC]+1
				direction[xk,yk]=kd[k,3]
				basins[xk,yk]=basins[xC,yC]
				count=count+1
				#nqueue=nqueue+1
			}
		}	
	#t4=proc.time()	
}
#print(paste(count, "available D4 neighbors found", nrow(queue1)))
#t0=proc.time()
#queue2=rbind(queue1, c(xk, yk, demtemp))
#t1=proc.time()
#t1-t0

##############
#Remove from the queue and move on
#t5=proc.time()
#nqueue=length(queue)/3

nqueue=length(queue1)/3
nqueuetemp=nqueue2+length(queuetemp)/3
#t6=proc.time()
#if you have 2 or more items in Q1 then you will have 1 left after
#deleting the current cell so continue with Q1
if(nqueue>1){
	#print("here")
	queue1=queue1[-pick,]
	nqueue=length(queue1)/3
	q1max=max(nqueue, q1max)
	#print(q1max)
#if you have only 1 item left in Q1 then it will be empty after deleting the
#current cell so see if you are don or if you need to grab another chunk from Q2
} else {
	#look and see if there are still values in Q2 to merge
	if(nqueuetemp>nchunk){
		t1=proc.time()
		split=split+1
		print(paste('Split:', split, 'Q2', nqueuetemp,  'nstep=', nstep, "Q1 Max:", q1max, "time", round((t1-t0)[3],2)))
		queue2=rbind(queue2,queuetemp)
		queuetemp=NULL
		qsort=queue2[order(queue2[,3]),]
		queue1=qsort[1:nchunk,]
		queue2=qsort[-(1:nchunk),]
		th=queue2[1,3]
		nqueue=length(queue1)/3
		nqueue2=length(queue2)/3
		q1max=0
		t0=proc.time()
		
	} else if(nqueuetemp<=nchunk & nqueue2 >0){
		split=split+1
		print(paste('Split:', split, 'Q2', nqueuetemp,  'taking last chunk, nstep=', nstep, "Q1 Max:", q1max))
		queue1=rbind(queue2, queuetemp)
		queuetemp=NULL
		queue2=NULL
		th=max(dem)*1.1
		nqueue=length(queue1)/3
		q1max=0
		nqueue2=0
	} else {
		print(paste('Q1 depleted, Q2', nqueue2,  'done!!'))
		nqueue=0
	}
}
#t7=proc.time()
#if there is only one row it will treat it as a vector and things will break
if(nqueue==1){
	queue1=matrix(queue1, ncol=3,byrow=T)
	}

nstep=nstep+1
if(printstep){print(paste("Step:", nstep, "NQueue:", nqueue, "Queue2:", nqueue2, "Qtemp:",length(queuetemp)/3 ))}

}



output_list=list("dem"=demnew, "mask"=mask, "marked"=marked, "step"= step, "direction"=direction, "basins"=basins)
return(output_list)

} # end function

