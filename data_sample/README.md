The files in this folder are sample data for the distortion analysis code.

The dataset includes the following files:

- `points.txt`: A text file containing pairs of control points, kindly provided by Dr. Martin Rickenbacher, martin.rickenbacher[at]bluewin.ch. The first column is the index, the second and third the old map coordinates and the fourth and fifth the new map coordinates.
- `old_map.jpg` and its 'world file' `old_map.jpgw`: An old map to analyze, named *"Die Landschaft Basel und das Frickthal - entworfen und mit beweglichen Typen gesetzt von W. Haas"*. It shows part of northern Switzerland and was produced in 1798 by W. Haas. The map has been scanned by swisstopo, Kartenarchiv (Signature LT K BL 1798).
![](./old_map.jpg?raw=true =250x)
- `new_map.gif` and its 'world file' `new_map.gfw`: A new map of the same area by swisstopo (scale 1:200,000 - projection CH1903 / LV03 EPSG::21781). Copyright 2005 swisstopo (JD052588).  Please note that part of it has been blanked for copyright reasons. You do not have the right to distribute, publish, or use the map for any other purposes.  
![](./new_map.gif?raw=true =250x)

As explained in the top-level [README](../README.md) file, these maps and data were kindly provided by other people and institutions. Please cite them accordingly.

This data is similar to the sample data set for [MapAnalyst](<http://mapanalyst.org/>) distortion analysis software, with one exception:
the coordinates of the old map were multiplied by 1000000, to prevent issues with rounding errors in the QGIS 'interpolate' routine.

 
