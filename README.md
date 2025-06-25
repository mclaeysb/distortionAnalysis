# distortionAnalysis

by Manuel Claeys Boùùaert

> [!NOTE]
> This distortion analysis approach, first published in 2015, has now, as of 2025, been implemented in the [Allmaps](https://allmaps.org/) project, offering a zero-setup, user friendly interface that greatly simplifies not only georeferencing maps but also inspecting their accuracy using differential distortion analysis.
> 
> Allmaps work with any map that is served as a [IIIF](https://iiif.io/) Resource. If your map is not available through IIIF yet, consider to use IIIF in your institution, or upload your map to a free IIIF server.
> 
> Use [Allmaps Editor](https://editor.allmaps.org/) to georeference maps and [Allmaps Viewer](https://viewer.allmaps.org/) to view the resulting [georeference annotation](https://preview.iiif.io/api/georef/extension/georef/). The Viewer application also allows for altering the transformation type (e.g. from *polynomial* to *thin plate spline*) and visualising distortions in the same way as in this codebase: using displacement vectors, distortion grid and *differential* approaches such as area and angular distortions. It currently only implements a *old-to-modern* comparison, showing the distortion of the old-to-modern transformation on the current, modern map.
>
> More information on [allmaps.org](https://allmaps.org)

## What it does

This code was built to analyse distortion in old maps.

Old, historical maps can be inaccurate for many reasons (primitive techniques, surveying errors, mixed inputs, conservation, etc). 
A distortion represents a loss in geographic information, but, by its specific origin, yields new historical information. 
Hence the study of distortions in old maps.

This code takes as input a set of common ground control points, identified both on the old map and on a modern reference map.
Based on these points, an interpolation is constructed between both maps. 
From this interpolation, the distortions are computed using common distortion visualisation techniques, including displacement vectors, distortion grids, Tissot indicatrices and a new technique called **Differential Distortion Analysis**.

The output consists of a set of `.WKT` files representing the displacement vectors, distortion grid and Tissot indicatrices, and a `.csv` file with the values of distortion measures, such as the **area distortion** `log2sigma` and **angular distortion** `2Omega`, computed at a number of mesh points and ground control points.

For further interpretation, visualise the output documents in a (web-)GIS application. For optimal visualisations, interpolate the point-wise results of the `.csv` file to create a continuous distortion raster image.
As an example: the visualised outputs for the sample data look like this: 
 
**Displacement vectors**  
![](/example_output/old_map_displacementVectors.png?raw=true)  
**Distortion grid**  
![](/example_output/old_map_distortionGrid.png?raw=true)  
**Indicatrices of Tissot**  
![](/example_output/old_map_indicatrices.png?raw=true)  
**Differential Distortion Analysis: area distortion `log2sigma`**  
![](/example_output/old_map_log2sigma_scalebar.png?raw=true)  
**Differential Distortion Analysis: angular distortion `2Omega`**  
![](/example_output/old_map_2Omega_scalebar.png?raw=true)  

Bonus: this code can also be used in similar settings where 2D objects are warped, such as abstracted maps (e.g. transit maps), mental maps, maps of deformed terrain (glaciers, tectonics), etc. Feel inspired!

## How to use it

The `.m` code is fully Matlab and Octave compatible.

In order to build the interpolation, you should have access to either the Octave splines package or Matlab CurveFitting Toolbox.
(The Octave splines package can be installed using `pkg install -forge -global -auto splines`)

The `example.m` file demonstrates how to load data, set optional parameters and call the main function `distortionAnalysis.m`.

## How to cite it

### Code and publication

Manuel Claeys Boùùaert wrote this code in 2015 for a project analysing distortion in the 1777 '[Carte de Cabinet](http://belgica.kbr.be/nl/coll/cp/cpFerraris_nl.html)' of the Austrian Netherlands by count de Ferraris.

If you want use this code for scientific work, **please refer to the following publication**:

*Claeys Boùùaert, M., De Baets, B., Vervust, S., Neutens, T., De Maeyer, P., Van De Weghe N., 2016. Computation and visualisation of the accuracy of old maps using differential distortion analysis. International Journal of Geographical Information Science.*  
<http://dx.doi.org/10.1080/13658816.2015.1127377>

The 'Author’s Original Manuscript (AOM)' of this article is also part of this repository.

If you want to learn more about the application of this technique on Farraris's 1777 Carte de Cabinet, read the following article:

*Vervust, S., Claeys Boùùaert; M., De Baets; B., Van de Weghe, N., De Maeyer, P., 2018. A Study of the Local Geometric Accuracy of Count de Ferraris's Carte de cabinet (1770s) Using Differential Distortion Analysis. The Cartographic Journal*  
<https://doi.org/10.1080/00087041.2017.1323148>

### Sample data

The sample data map is courtesy of swisstopo (map collection, LT K BL 1798).   
The sample data ground control points are courtesy of Martin Rickenbacher (Bern).

This sample data is similar to the one used in the [MapAnalyst](<http://www.mapanalyst.org>) software by Bernhard Jenny, which has inspired us during our work.

### Octave splines package

The [splines package](http://octave.sourceforge.net/splines) for GNU Octave is written by Nir Y. Krakauer and others. It features the commands used in this code since version 1.2.9.