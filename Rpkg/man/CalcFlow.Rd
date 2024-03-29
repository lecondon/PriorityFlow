% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CalcFlow.R
\name{CalcFlow}
\alias{CalcFlow}
\title{Calculate Overland Flow}
\usage{
CalcFlow(
  file.path,
  run.name,
  file.nums,
  slopex.file,
  slopey.file,
  overland.type,
  mannings,
  epsilon = 1e-05,
  dx = 1,
  dy = 1,
  mask
)
}
\arguments{
\item{file.path}{path where pressure and slope files are located, this is also where flow files will be written to}

\item{run.name}{ParFlow run name this is used to read the pressure files assumign they are named with the
standard ParFlow naming convention (i.e. runname.out.press.00000.pfb)}

\item{file.nums}{list of file numbers to be processed}

\item{slopex.file}{name of the x slopes file  (note this must be a pfb file)}

\item{slopey.file}{name of the y slopes file  (note this must be a pfb file)}

\item{overland.type}{type of overland flow used for ParFlow simulation. Choices are 
"OverlandFlow" or "OverlandKinematic" this shoudl match what you used in your tcl script.}

\item{mannings}{manning roughness coefficient, can be either a value or a matrix}

\item{epsilon}{Epsilon used only in the OverlandKinematic formulation, this is set with 
the Solver.OverlandKinematic.Epsilon key in ParFlow (default value 1E-5)}

\item{dx}{grid size in x direciton [l] (default value 1)}

\item{dy}{grid size in y direciton [l] (default value 1)}

\item{mask}{mask with ones for cells to be processed and zeros for everything else - defaults to a mask of all 1's
if a mask is a provide it shoudl be a [nx,ny] matrix such that such that [1,1] is the lower left corner of the 
domain and [nx,ny] is the upper right}
}
\description{
Function to calculate overland flow from ParFlow pressure file outputs. This function will write three pfb outputs 
for every timestep (1) outflow - this is the volumetric outflow from each grid cell [l3/t], (2) q_east - the 
volumetric flow across the east face of each cell (Note: this will have dimensions nx+1, ny to account for the western edge
of the domain) (2) q_north - the volumetric flow across the north face of each cell (Note this will have dimensions nx, ny+1
 to account for the southern edge of the domain)
}
