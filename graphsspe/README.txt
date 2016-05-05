PRINCIPLE
================

The fishing arena is discrete and 
is defined as a graph (a set of vertices and edges)
to facilitate the modelling of individual vessel 
movements and directed individual decisions 
by enabling the use of well-known 
data-tree-structure algorithms. 
In particular, Dijkstra's (1959) algorithm 
is used to search for the shortest path(s)
between a departure and a destination area 
(in our case, from a port to a fishing ground and vice versa). 
Nevertheless, the movements of vessels travelling 
along the graph are still continuous, 
i.e., vessels are not necessarily stopping on each node,
 to keep unbiased travelling distance estimates. 
The graph nodes lie on a regular grid, and each node 
is linked with its eight respective neighbour nodes 
to enable vessel movement in eight geographical directions 
(North, N-E, East, S-E, South, S-W, West, N-W). 

For example, for the 'balticonly' application,
the entire graph has three levels of grid resolution 
with the ruling principle that one management area can only 
have one grid resolution to avoid bias for underlying stocks 
that are by definition given by the management area.
The graph nodes are spaced accordingly by 10 km, 25 km, and 80 km
from the finest to the broader resolution scale (grid) 
for the Skagerrak/Kattegat and Western Baltic areas, 
for the North Sea and the Central Baltic Sea areas, 
and for the areas outside, respectively. 

R PSEUDOCODE
=================

The graph is generated with an R routine:

1- define a regular grid of points (nodes) over decimal latitude on longitude
   (**TO DO: if the area is too extensive, risk of distance bias**, then
   need to change for a grid of equi-distant points e.g. based on UTM coordinates)

2- clip points with polygons, usually with polygons from a GIS shape file
   (for example all the polygons land to only keep the points at sea)

2bis- As a refinment, a small dilatation of the coastline can be apply 
       to avoid to close points on the shore

3- Harbour nodes are added

4- find the 8 neighbours to each node (the R function 'knearneigh' in 'spdep' R package)
   nbk           <- 8
   graph.neightb <- knearneigh(coord[,c("x","y")], k=nbk, longlat = TRUE)
   dd            <- nbdists(knn2nb(graph.neightb) , coord[,c("x","y")], longlat = TRUE)  # return distance in km
   node          <- rep(1:nrow(graph.neightb$nn), each=nbk)
   neighbour     <- matrix(t(graph.neightb$nn), ncol=1)
   dist.km       <- unlist(dd)
   graph         <- cbind(node, neighbour, dist.km)

5- reduce number of neighbours specifically for the harbour nodes
 in order to avoid potential weird connections for example across land

6- format and save for suitable use in c++
i.e. coordXX.dat storing the x (i.e. long), y (i.e. lat) coord of the node 0 to N-1 and 
graphXX.dat storing the edges (from,to) and the distance in km between from,to
(everything in one column, so you also need the number of nodes and edges)



DIFFERENT TYPES OF GRAPHS
=================

Child graphs, (if parent is 11, then children are 111, 112, 113, etc.)
can be further derived by delimiting polygons assigning a penalty distance
to some of the edges in order to avoid the pathfinding algo to pass through them
anymore. Those polygons are also usually obtained from GIS shape files.


ADDITIONAL MINOR FILES
==================

code_area_for_graphXX_points.dat
with an assigned code from a GIS layer, e.g. ICES areas
However, used in the core c++ model.

code_square_for_graphXX_points.dat
A more readable format of the coord, aslo with ICES rectangle coding
Not used in the c++ core model, just given for info
Furthermore the number of digits is less precise.

nodes_in_polygons_a_graphXX_quarterX.dat
Given for info,
used in the core c++ model.
Per quarter of the year






