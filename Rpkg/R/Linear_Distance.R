#' Calculate the minimum linear distance from a feature set
#' 
#' Calculates the minimum linear distance between every point in a 2D array
#' and a mask of target points. Note that this function assumes dx=dy if this is not the case it will not work. 
#' @inheritParams D4TraverseB
#' @param  target_points The matrix with a value of 1 for all cells that you would like to calculate the distance to and 0 for all other cells
#' @param  cell_size The size of a grid cell (**NOTE that this function assumes dx=dy, i.e. square grid cells)
#' @export

LinDist=function(target_points, mask, cell_size=1, printflag=F){

nx=nrow(target_points)
ny=ncol(target_points)

if(missing(mask)){mask=matrix(1, nrow=nx, ncol=ny)} #default to processing everything

# incomplete matrix for cells still needing to be processed (1= cells needing to be processed, 0= done cells)
# Target points and points outside the mask are initiliazed to 0
incomplete = (1-target_points)*mask

n_missing = sum(incomplete)
distance=matrix(0, nrow=nx, ncol=ny)
s=1

while (n_missing>0 & s<max(nx,ny)){
  print(paste("Step=", s, "N_missing=", n_missing))
  for(r in 0:s){
      d_temp = ((s*cell_size)^2+(r*cell_size)^2)^(0.5)
      print(paste("Rotation=", r, "N_missing=", n_missing, "distance=", d_temp))
      #counting up the number of target points within a given step (s= later distance from center) and rotation (r=distance off vertical or horizonatl axes)
      temp_count = matrix(0, nrow=nx, ncol=ny)
      if(s<nx & r<ny){
        temp_count[1:(nx-s), 1:(ny-r)] = temp_count[1:(nx-s), 1:(ny-r)] + target_points[(s+1):nx, (r+1):ny ]
        temp_count[1:(nx-s), (r+1):ny] = temp_count[1:(nx-s), (r+1):ny] + target_points[(s+1):nx, 1:(ny-r)]
        temp_count[(s+1):nx, 1:(ny-r)] = temp_count[(s+1):nx, 1:(ny-r)] + target_points[1:(nx-s), (r+1):ny ]
        temp_count[(s+1):nx, (r+1):ny] = temp_count[(s+1):nx, (r+1):ny]  + target_points[1:(nx-s), 1:(ny-r)]
      }
      if(r<nx & s<ny){
        temp_count[1:(nx-r), 1:(ny-s)] = temp_count[1:(nx-r), 1:(ny-s)] + target_points[(r+1):nx, (s+1):ny ]
        temp_count[1:(nx-r), (s+1):ny] = temp_count[1:(nx-r), (s+1):ny] + target_points[(r+1):nx, 1:(ny-s)]
        temp_count[(r+1):nx, 1:(ny-s)] = temp_count[(r+1):nx, 1:(ny-s)] + target_points[1:(nx-r), (s+1):ny ]
        temp_count[(r+1):nx, (s+1):ny] = temp_count[(r+1):nx, (s+1):ny] + target_points[1:(nx-r), 1:(ny-s)]
      }
      #image.plot(temp_count)
      temp_count[temp_count>0] = 1
      # record the distance for any cell with a temp count=1 that is still incomplete
      distance = distance + temp_count * incomplete * d_temp
      
      #set all the cells that had a temp_count=1 to complete
      incomplete= incomplete-temp_count
      incomplete[incomplete<0]=0
      
      #count up the number of cells still needing a distance
       n_missing = sum(incomplete)
       
       #stop if you ar out of missing cells
       if(n_missing==0){break}
  } #end for rotation
  s = s+1
} #end while
distance[mask==0]=NA

return(distance)
}
