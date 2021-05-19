# flow_direction_generalization
Python code implementing COTAT<sup>1</sup>, NSA<sup>2</sup> and DMM<sup>3</sup> methods for generalization of flow direction rasters on continental to global scale

NSA and DMM require flow accumulation raster derived from flow directions. COTAT method uses both flow direction and flow accumulation rasters.
All three methods have a parameter of cell factor: integer cellsize multiplication factor ***k***, i.e. the number by which input cellsize is multiplied in output raster. Cell areas in output model are ***k<sup>2</sup>*** times larger. COTAT has an additional parameter: ***∆A<sub>thr</sub>*** in **km<sup>2</sup>** which defines when the tracing on input high resolution flow direction graph is stopped and coarse resolution flow direction is assigned. Empirical equation for ***∆A<sub>thr</sub>*** in **km<sup>2</sup>** based on output cellsize ***L*** in degrees: ***∆A<sub>thr</sub>=22500L<sup>2</sup>***. Classes COTAT, NSA and DMM that run corresponging methods are imported from methods.py

In order for parameter ***∆A<sub>thr</sub>*** to make sense when dealing with rasters in geographic coordinate system, flow accumulation for COTAT method must be calculated using cell area raster as weight raster in flow accumulation tool (ArcGIS, GRASS etc). Cell area raster is calculated as spheroidal trapezoid areas using **CellAreas** class from **calculate_cell_areas.py** taking sample raster as an input (initial high resolution flow direction raster can be passed here).

Metrics are calculated for each watershed separately. Metrics include watershed area (**km<sup>2</sup>**), orthogonal to diagonal flow length ratio, and mean flow length. Metrics for watersheds are calculated using **WatershedStats** class from metrics.py. Flow direction raster and watershed raster derived from flow directions are required to calculate metrics. Results are saved to a **.csv** file.

Dependencies: gdal, numpy, heapq, geopy. 

1. Reed S.M. Deriving flow directions for coarse-resolution (1-4 km) gridded hydrologic modeling // Water Resources Research. 2003, vol. 39, iss. 9, pp. 1–11.
2. Fekete B.M., Vörösmarty C.J., Lammers R.B. Scaling gridded river networks for macroscale hydrology: Development, analysis, and control of error // Water Resources Research. 2001, vol. 37, iss. 7, pp. 1955–1967.
3. Olivera F. et al. Extracting low-resolution river networks from high-resolution digital elevation models // Water Resources Research. 2002, vol. 38, iss. 11, pp. 1-8.
