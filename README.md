# topographic-processing
This repo contains a set of R Functions and example workflow scripts for processing DEMs to create input for hydrologic models

### DEM processing
The DEM processing code is a modified version of the 'Priority Flood' algorithm which is a depression filling algorithm (Wang and Liu 2006; Barnes et al., 2014; Zhou, 2016).  This is an optimal approach that processes the DEM to ensure every cell drains to the border.

As implemented here there are options to ensure drainage to the edges of a regular rectangular domain or the user can provide a mask for an irregular domain.

Additionally, a second processing option is provided if there is an a-priori river network that you would like to enforce. In this case the river network is provided first as a mask to the processing algorithm to ensure that every identified river cell drains to the boundary of the domain (regular or irregular) next the remaining cells are processed using the river network as the 'boundary', ensuring that every other cell drains either to a river or to a boundary. For examples of this approach refer to 'Workflow_Example3.R' and 'Workflow_Example4.R'

### Slope calculations
The slope function 'Slope_Calc_Upwind.R' calculates slopes in the x and y direction up-winded to be consistent with ParFlow. The most basic application of the function will calculate slopes in the x and y directions for the entire domain based on the input DEM

For more advanced processing there are the following options (refer to the function for details on how to implement flags):
+ Maximum slope threshold: A maximum absolute value slope can be set
+ Scaling secondary flow directions: By default slopes are calculated in the x and y direction. Using the flow direction file the slope perpendicular to the primary flow direction can be scaled to some fraction of the primary flow direction. If the fraction is set to 0 this will result in slope files with only slopes in x or y for a given cell (similar to previous GRASS processing approaches)
+ Enforcing primary flow directions along river cells: The flow scaling applied above occurs over the entire domain. In addition, if a mask of river cells is applied one scaling value can be applied to the rest of the domain while flow in the primary directions is enforced along the specified river cells (i.e. scaling ratio of 0 for river cells).
+ Border cells: By default all border cells are processed to point out. An optional border direction mask can be input to set some borders to point in (note where overland flow occurs rivers will still be pointed out)

### Workflow scripts
Four workflow scripts are provided to demonstrate processing for four examples of increasing complexity
1. 'Workflow_Example1.R': Rectangular domain with no river network
2. 'Workflow_Example2.R': Irregular domain with no river network
3. 'Workflow_Example3.R': Rectangular domain with river network
4. 'Workflow_Example4.R': Irregular domain with river network
