PriorityFlow
=======
*PriorityFlow* is a toolkit for topographic processing for hydrologic models. This repo contains an R package and a set of workflow examples (see instructions below).  

#### Development Team
+ Laura Condon (lecondon@email.arizona.edu)
+ Reed Maxwell (rmaxwell@mines.edu)

#### Citation
For more details on the model and if you use PriorityFlow in published work please cite the following reference:  
   *Condon, LE and RM Maxwell (2019). Modified priority flood and global slope enforcement algorithm for topographic processing in physically based hydrologic modeling applications. Computers & Geosciences, [doi:10.1016/j.cageo.2019.01.20](https://doi.org/10.1016/j.cageo.2019.01.020).*

Installing the R-package
--------------------
The priority flow functions are provided as an R package. To install it you will need to do the following:
1. To use R packages from GitHub you need to have the devtools package installed:
```
install.packages('devtools’)
library(devtools)
```

2. Next you can install the PriorityFlow package:
```
install_github("lecondon/PriorityFlow", subdir="Rpkg”)
library('PriorityFlow’)
```
3. If its installed correctly you should be able to run help and get information. Clicking on the index from here will give you a list of all of the available functions.
```
help('PriorityFlow’)
```
*Note: If you are using the R-package then you can get rid of the 'source' calls at the top of the example workflows because the functions will already be loaded with the library.

Getting  started
--------------------
The best way to get started with this toolset would be to walk through the 'Workflow_Example' R notebook or html file in this directory. Additionally, you can refer to the help files for all of the functions individually in index when if you run 'help('PriorityFlow')'

DEM Processing
--------------------
The DEM processing code is a modified version of the 'Priority Flood' algorithm which is a depression filling algorithm (Wang and Liu 2006; Barnes et al., 2014; Zhou, 2016).  This is an optimal approach that processes the DEM to ensure every cell drains to the border.

As implemented here there are options to ensure drainage to the edges of a regular rectangular domain or the user can provide a mask for an irregular domain. NOTE: if you are providing a mask for an irregular domain boundary you must ensure that your mask is D4 contiguous. See the mask tip below for one approach to achieve this using grass.

Additionally, a second processing option is provided if there is an river network that you would like to enforce. In this case, the river network is provided first as a mask to the processing algorithm to ensure that every identified river cell drains to the boundary of the domain (regular or irregular). This step will also ensure a D4 connected river network (i.e. stair stepping around any diagonal river neighbors). Next the remaining cells are processed using the river network as the boundary, ensuring that every other cell drains either to a river or to a boundary. For examples of this approach refer to `Workflow_Example3.R` and `Workflow_Example4.R`.

Slope Calculations
--------------------
There are two slope calculations functions in this repo:

`Slope_Calc_Upwind.R`  calculates slopes in the x and y direction down-winded to be consistent with the ParFlow OverlandFLow boundary condition.

`Slope_Calc_Standard.R` calculates slopes in the x and y direction using indexing to be consistent with the ParFlow OverlandKinematic and OverlandDiffusive boundary conditions. This is the approach that is used in the main workflow example.


Workflow Scripts
--------------------
1. `Workflow_Example.Rmd`: This is the most updated  workflow example and  the  one  I recommend starting from.

The next four examples show the older slope calculation function with downwinding
1. `Downwinding_Workflow_Example1.R`: Rectangular domain with no river network
2. `Downwinding_Workflow_Example2.R`: Irregular domain with no river network
3. `Downwinding_Workflow_Example3.R`: Rectangular domain with river network
4. `Downwinding_Workflow_Example4.R`: Irregular domain with river network


Tips
--------------------
If you want to process your DEM within a pre-defined watershed mask and you need help creating that mask. An example workflow using QGIS and GRASS to ensure a D4 connected mask (i.e. one where you don't have any cells that are only connected to the rest of the domain diagonally):
1. create the mask in QGIS
2. Using GRASS, clump it: r.clump
3. change no data value to 0: r.null
4. identify the ID of the big contiguous block (clump) you want to keep
5. re-mask:`("NullRaster@1"=65)*1+("NullRaster@1" != 65)*0`

License
--------------------
Copyright (C) 2018  Laura Condon

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation version 3 of the License

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>

References
--------------------
+ Barnes, R., C. Lehman, and D. Mulla, Priority-flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models. Computers & Geosciences, 2014. 62: p. 117-127.
+ Wang, L. and H. Liu, An efficient method for identifying and filling surface depressions in digital elevation models for hydrologic analysis and modelling. International Journal of Geographical Information Science, 2006. 20(2): p. 193-213.
+ Zhou, G., Z. Sun, and S. Fu, An efficient variant of the Priority-Flood algorithm for filling depressions in raster digital elevation models. Computers & Geosciences, 2016. 90: p. 87-96.
