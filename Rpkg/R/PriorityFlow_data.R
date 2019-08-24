#' Example watershed DEM for PrioirtyFlow 
#' 
#' This is a small elevation dataset (215km by 172 km at 1km spatial resolution) 
#' 
#' @docType data
#' @usage data(DEM)
#' @format a matrix of elevaiton values with nrow=nx and ncol=ny dimensions of the domain
#' @keywords datasets
"DEM"

#' Watershed Mask
#'
#' A mask showing the watershed drainage area for the test domain
#' 
#' @format a matrix of 0's and 1' showing the watershed extent (1=inside the watershed, 0=outside the watershed), 
#' where nrow=nx and ncol=ny for the domain.
#' @keywords datasets
"watershed.mask"

#' River Mask
#'
#' A mask showing and example river network for  the test domain
#' 
#' @format a matrix of 0's and 1' showing the location of the river cells (1=river, 0=nonriver), where nrow=nx
#' and ncol=ny for the domain.
#' @keywords datasets
"river.mask"