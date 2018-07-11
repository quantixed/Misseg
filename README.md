# Misseg

A set of tools to analyse the spatial organisation of features in a cell and also measure their intensity.

In FIJI, the Cell Counter is used to highlight specific features of interest (and select some background regions). XML files from Cell Counter are saved to a directory.

TIFF files of all the images are created - if necessary using `dvExtractor.ijm` and saved to a directory.

In Igor, either the spatial or intensity analysis is done using the appropriate menu item.

## Spatial analysis

For spatial analysis, point Igor at your directory of XML files and it will do the rest. Currently this is hard-coded to our specific use case.

## Intensity analysis

For intensity analysis a simple panel allows the user to select the XML and TIFF directories and then set other parameters. These include what is imaged in each channel, which objects have been picked with the cell counter, the voxel size and the number of channels. An important setting is the "radius" which determines the size of the sphere used to analyse the mean intensity centred on the points picked out.

This analysis in a spherical ROI is similar to the function in IMARIS. It works using extraction of a hypercube of data and multiplication of a binarised spherical cube.