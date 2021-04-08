# spain_popgen

Feel free to use these scripts for academic or personal use.  These are based on work for Clare Bycroft's PhD thesis, and subsequent [publication](https://doi.org/10.1038/s41467-018-08272-w) of a study in Spanish genetic population structure. 

DISCLAIMER: Scripts may contain hard-coded directories, or will try to load files that you don't have access to. I recommend copy-pasting sections of scripts where required for your own code.

*Please cite:*

Bycroft, C., Fernandez-Rozadilla, C., Ruiz-Ponte, C. et al. Patterns of genetic differentiation and the footprints of historical migrations in the Iberian Peninsula. Nat Commun 10, 551 (2019). https://doi.org/10.1038/s41467-018-08272-w


## How to make pretty maps 
The following outlines how to build similar maps to Figures 1 and 4 in the [publication](https://doi.org/10.1038/s41467-018-08272-w)) - for any arbitrary map, and any set of sampled point locations.

### Function library

I have a large function library, [Other/ClaresFunctions.R](https://github.com/cgbycroft/spain_popgen/blob/master/Other/ClaresFunctions.R), which is basically my PhD workhorse, developed over 4 years ;D.  If other scripts don't run, it's likely because they're calling a function from ClaresFunctions.R, so if you can get it so you can source the library script it might make other things easier, as you won't have to dig around the large library so much.

You won't be able to directly run `source("Other/ClaresFunctions.R")` as it is, because it reads in some objects from local files you won't have access to.  However, I think most of the lines you'll have to remove that load in data are at the beginning on line 85, but keep lines 93-96 and 99 as these load in maps that I've included in the repository, and which are used further down in the script.


### Mapping fineSTRUCTURE clusters

I have included a short tutorial for mapping using the spplot package, which you might find helpful (uses a map of New Zealand): [mapping/genericMapping.R](https://github.com/cgbycroft/spain_popgen/blob/master/mapping/genericMapping.R)

More specifically, the mapping of fineSTRUCTURE clusters is done with the two scripts (the first one you really only want to run once):

[densityMaps/getSpatialGrid_v2.R](https://github.com/cgbycroft/spain_popgen/blob/master/densityMaps/getSpatialGrid_v2.R), which generates a spatial grid and compute distances between samples, and then

[chromopainter/GaussianMappingGeneric.R](https://github.com/cgbycroft/spain_popgen/blob/master/densityMaps/continuousGaussianPlotting.R)

which generates the layers for plots.

In both script you'll have to replace the object "espMapTotalProj" with your own SpatialPolygonsDataFrame object, which contains polygons of the map you want to plot.  You will also need the X and Y coordinates (in longitude/latitude) of the samples you want to plot (you'll have to adapt lines 69 - 84 of getSpatialGrid_v2.R).   I've added R objects of the maps I used to the folder "mapping", so if you want to play around you can load in those (lines 1473-1586 in Other/ClaresFunctions.R creates the projected (\*Proj) map and some other mapping objects).

For the GaussianMappingGeneric.R script you'll need the finestructure results as a data frame called "split.matrix", with each row an individual, and each column a set of cluster assignments based on successive layers of the tree.  You can form this matrix from the finestructure results in a few ways. I've done it by just slicing the tree at each branch split keeping the information in the lengths of the branches.  I've written two functions to do this based on the XML files.  See Other/ClaresFunctions.R

line 3547:  `readFine()`

and

line 3595:  `newSplitMatrix()`

These functions call some other functions, so you might have to dig them out of the same library script (they should be nearby), or just source("Other/ClaresFunctions.R") if you can. 

The function `plotPointsOnMaps2()` is useful for plotting coloured points on a map.  `plotDenisty()` and `plotDenistyContinuous()` are used for plotting density maps (note the typo that I never got around to fixing). 



### Spatially smoothed ancestry profiles

Once you have your map object, a grid object, your sample locations, a physical distance matrix ([densityMaps/getSpatialGrid_v2.R](https://github.com/cgbycroft/spain_popgen/blob/master/densityMaps/getSpatialGrid_v2.R)), and you have a matrix of coancestry vectors, you can compute the spatially smoothed ancestry profiles. 

You will next need to get a matrix of weights for each grid-point/sample pair ([densityMaps/getVariableBandwidth.R](https://github.com/cgbycroft/spain_popgen/blob/master/densityMaps/getVariableBandwidth.R)).

Then the spatial ancestry profiles are computed in the script mixtureModel/Plot_mixture_model2.R.  You can basically ignore everything until line 188-253. You'll need two other objects: "SUMMARYmat_1" which a nxK matrix of coancestry vectors for each recipient individual (chunklengths from each donor group summed), and SUMMARYmat is a KxK matrix with the coancestry vectors for each donor group (averaged within each donor group).

I created these objects using the function `getCountSummaries()` in [Other/ClaresFunctions.R](https://github.com/cgbycroft/spain_popgen/blob/master/Other/ClaresFunctions.R).

You can plot the output using the function `plotDenistyContinuous()`.  For an example of how to plot it with a map and also sample locations see lines 298 - 367 of Plot_mixture_model2.R.




### Plotting spatially smoothed continuous variables (generic)

The script continuousGaussianPlotting.R does plotting of any continuous variable as a smoothed gradient on a map.  However, if you have already run 'getVariableBandwidth.R' it is actually a bit redundant, as you just need to then run the following lines:

```
sumd <- colSums(d)

out = inputVector*d
toPlot <- rowSums(t(out)/sumd)
plots = plotDenistyContinuous(toPlot)
```

where `d` = distance-based contributions of samples to each grid point (n-samples x number-of-grid-points) calculated in getVariableBandwidth.R, or by any other means (see continuousGaussianPlotting.R for other possible examples); and inputVector = a vector of values for some continuous variable (one element per sample), in the order of the rows of d.




