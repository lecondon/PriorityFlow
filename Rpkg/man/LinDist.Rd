% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Linear_Distance.R
\name{LinDist}
\alias{LinDist}
\title{Calculate the minimum linear distance from a feature set}
\usage{
LinDist(target_points, mask, cell_size = 1, max_dist = Inf, printflag = F)
}
\arguments{
\item{target_points}{The matrix with a value of 1 for all cells that you would like to calculate the distance to and 0 for all other cells}

\item{mask}{Mask with ones for cells to be processed and zeros for everything else - defaults to a mask of all 1's}

\item{cell_size}{The size of a grid cell (**NOTE that this function assumes dx=dy, i.e. square grid cells)}

\item{max_dist}{Maximum distance to check. All cells with no distance found <= max_dist will be assigned max_dist}

\item{printflag}{optional flag to print function progress. this is a test}
}
\description{
Calculates the minimum linear distance between every point in a 2D array
and a mask of target points. Note that this function assumes dx=dy if this is not the case it will not work.
}
