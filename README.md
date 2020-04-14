# Misseg

A set of tools to analyse the spatial organisation of features in a cell, measure their intensity and generate averaged views of these features.

In FIJI, the Cell Counter is used to highlight specific features of interest (and select some background regions). XML files from Cell Counter are saved to a directory.

TIFF files of all the images are in a directory to be loaded in if required (they can be created if necessary using `dvExtractor.ijm`).

In Igor, an analysis type is done using the appropriate menu item.

## Spatial analysis

### Standard spatial analysis

For standard spatial analysis, point Igor at your directory of XML files and it will do the rest. No Tiffs are needed. The points will be scaled to the position of the two spindle poles so that a normalised view of all points in the dataset can be visualised and compared.

### Exclusion zone spatial analysis

This workflow shows the distribution of pointes relative to the exclusion zone. It needs a set of segmented images to be provided. Depending on the data, these can be automatically generated or if necessary manually segmented. The names must match the Cell Counter files.

### Hybrid method

We are using a hybrid between the two analysis types above. Select the checkbox from the dialog that shows after the Exclusion Zone method is selected.

## Image analysis

### Intensity analysis

For intensity analysis, a simple panel allows the user to select the XML and TIFF directories and then set other parameters. These include what is imaged in each channel, which objects have been picked with the cell counter, the voxel size and the number of channels. An important setting is the "radius" which determines the size of the sphere used to analyse the mean intensity centred on the points picked out.

This analysis in a spherical ROI is similar to the function in IMARIS. It works using extraction of a hypercube of data and multiplication of a binarised spherical cube.

### Image clipping

Image clipping will excise a little clip of image data centred on each specified point. A simple panel allows the user to select the XML and TIFF directories and then set other parameters.

These clips are rotated so that the spindle pole (object 2) is directly below the object. The clips are padded so that the contour point stays in the centre of the clip. These rotated clips are averaged. A subsampling method is used to make the averages from an equally sized sample.

Note that the furthest half of object 5 is used as the pool for averaging.

--

Note that flexibility of these analyses is limited. Further modifications to the code are needed to (for example) analyse more or less than 4 object classes.

## Kinetochore coloc

A series of scripts to examine co-localisation of kinetochore signals in 3D image stacks.

- `kinetochoreColoc.ijm` loads dv files segments them and uses 3DOC to find the positions of two signals in channels 3 and 4. This is a substitute for Imaris (since we were having problems with the bridge to MATLAB and we did not have a current version)
- `kinetochoreSegmentation.ijm` is a tool to check how the segmentation (thresholding) performed in the previous script (it is not needed for analysis _per se_)
- `dv2RGB.ijm` will convert a folder of dv files to RGB tiffs for classification in Igor
- `kinetochoreColoc.ipf` is a menu-based workflow to classify the kinetochores from `kinetochoreColoc.ijm` and then analyse their proximity. There is an option to verify the classification done by a user and even to add additional kinetochores that have been missed by the automated procedure.
