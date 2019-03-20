PriorityFlow
=======
*PriorityFlow* is a toolkit for topographic processing for hydrologic models. This repo contains a set of R Functions and example workflow scripts.

#### Development Team
+ Laura Condon (lecondon@email.arizona.edu)
+ Reed Maxwell (rmaxwell@mines.edu)

#### Citation
For more details on the model and if you use EcoSLIM in published work please cite the following reference:  
   *Condon, LE and RM Maxwell (2019). Modified priority flood and global slope enforcement algorithm for topographic processing in physically based hydrologic modeling applications. Computers & Geosciences, [doi:10.1016/j.cageo.2019.01.20](https://doi.org/10.1016/j.cageo.2019.01.020).*


DEM Processing
--------------------
The DEM processing code is a modified version of the 'Priority Flood' algorithm which is a depression filling algorithm (Wang and Liu 2006; Barnes et al., 2014; Zhou, 2016).  This is an optimal approach that processes the DEM to ensure every cell drains to the border.

As implemented here there are options to ensure drainage to the edges of a regular rectangular domain or the user can provide a mask for an irregular domain. NOTE: if you are providing a mask for an irregular domain boundary you must ensure that your mask is D4 contiguous. See the mask tip below for one approach to achieve this using grass.

Additionally, a second processing option is provided if there is an river network that you would like to enforce. In this case, the river network is provided first as a mask to the processing algorithm to ensure that every identified river cell drains to the boundary of the domain (regular or irregular). This step will also ensure a D4 connected river network (i.e. stair stepping around any diagonal river neighbors). Next the remaining cells are processed using the river network as the boundary, ensuring that every other cell drains either to a river or to a boundary. For examples of this approach refer to `Workflow_Example3.R` and `Workflow_Example4.R`.

Slope Calculations
--------------------
The slope function `Slope_Calc_Upwind.R` calculates slopes in the x and y direction up-winded to be consistent with ParFlow. The most basic application of the function will calculate slopes in the x and y directions for the entire domain based on the input DEM

For more advanced processing there are the following options (refer to the function for details on how to implement flags):
+ Maximum slope threshold: A maximum absolute value slope can be set.
+ Scaling secondary flow directions: By default slopes are calculated in the x and y direction. Using the flow direction file which is produced in the DEM processing step, the slope perpendicular to the primary flow direction can be scaled so that they do not exceed some fraction of the primary flow direction. If the fraction is set to 0 this will result in slope files with only slopes in x or y for a given cell (similar to previous GRASS processing approaches).
+ Enforcing primary flow directions along river cells: The flow scaling applied above occurs over the entire domain. In addition, if a mask of river cells is applied one scaling value can be applied to the rest of the domain while flow in the primary directions is enforced along the specified river cells (i.e. applying a scaling ratio of 0 for river cells).
+ Smoothing slopes along river reaches: River reaches are calculated using the sub-basin function and can vary in size baed on the the area threshold set for sub-basins. Smoothing can be applied along river reaches either using the average slope for the cells in the river reach or the average slope of the entire subbasin.
+ Border cells: By default all border cells are processed to point out. This can be changed to have all border cells to point by changing *borderdir* option in the slope function to 2. In this case, border cells will point in, except wherever there is one or more cells draining to a border cell, this will pointed out to allow surface water bodies to exit. For more complicated domains, an optional border direction mask can be input to set some borders to point out and some in (Note that as with the previous setting  wherever there are one or more cells draining to a border cell the cell will pointed out to allow surface water bodies to exit regardless of what is specified in the mask).

Workflow Scripts
--------------------
Four workflow scripts are provided to demonstrate processing for four examples of increasing complexity
1. `Workflow_Example1.R`: Rectangular domain with no river network
2. `Workflow_Example2.R`: Irregular domain with no river network
3. `Workflow_Example3.R`: Rectangular domain with river network
4. `Workflow_Example4.R`: Irregular domain with river network


Tips
--------------------
An example workflow using QGIS and GRASS to ensure a D4 connected mask (i.e. one where you don't have any cells that are only connected to the rest of the domain diagonally):
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
