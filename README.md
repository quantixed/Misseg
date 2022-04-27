# Misseg

Code written for `ImageJ`, `R` and `IGOR Pro` to analyse data associated with manuscript `#p055p032`

**Endomembranes promote chromosome missegregation by ensheathing misaligned chromosomes**

Nuria Ferrandiz, Laura Downie, Georgina P. Starling, Stephen J Royle

*J Cell Biol*, (2022) 221: e202203021 doi: [10.1083/jcb.202203021](https://doi.org/10.1083/jcb.202203021)

--

## Misseg

The main program in this repo is `Misseg/Misseg.ipf`.
Full details of how to use it with `IGOR Pro` can be found in `Misseg/README.md`

The other directories are:

- [ERRemoval](#er-removal)
- [kinetochoreColoc](#kinetochore-coloc)
- [Progression](#progression)
- [SankeyFate](#sankey-fate)

## ER Removal

`CongressionAnalysis` - Quantifying the congression of a misaligned chromosome to the plate.

DNA channel is processed and segmented in ImageJ manually.
User defines the segments belonging to chromosome and plate using Magic Wand tool and adding to ROI manager.
ROIs are converted to csv using a python script and files are saved into a folder for chromosome and a folder for plate, within a folder for the movie.

`ER Removal` - Quantifying the clearance of ER in mitotic cells.

Using the same ImageJ scripts, the ER channel was processed and segmented using `ER_Removal.ijm` or the two subscripts (divided to allow the slowest step to be completed and different segmentation approaches tried with the other).
Resulting data are processed in Igor using `MissegERRemoval.ipf`


## Kinetochore coloc

A series of scripts to examine co-localisation of kinetochore signals in 3D image stacks.

- `kinetochoreColoc.ijm` loads dv files segments them and uses 3DOC to find the positions of two signals in channels 3 and 4. This is a substitute for Imaris (since we were having problems with the bridge to MATLAB and we did not have a current version)
- `kinetochoreSegmentation.ijm` is a tool to check how the segmentation (thresholding) performed in the previous script (it is not needed for analysis _per se_)
- `dv2RGB.ijm` will convert a folder of dv files to RGB tiffs for classification in Igor
- `kinetochoreColoc.ipf` is a menu-based workflow to classify the kinetochores from `kinetochoreColoc.ijm` and then analyse their proximity. There is an option to verify the classification done by a user and even to add additional kinetochores that have been missed by the automated procedure.

![img](img/coloc.png?raw=true "image")

The classifier above allows a user to scroll through a Z-stack and classify kinetochores as Aligned, Misaligned (called Free in the paper) or Misaligned - ER (called Ensheathed in the paper).

## Mad2

These experiments were analysed in Igor using `Mad2/ChromosomeCongression.ipf` by processing csv outputs generated in Fiji using `Mad2/ChromosomeCongressionMagicWand.ijm`

The method for chromosome and plate objects segmentation is as described for ER Clearance above.

## MnQuant

ImageJ scripts were used to generate the orthogonal views shown in the figure as well as do the quantification of micronuclei.
Outputs from the quantification scripts are in `Data/` for LBR and BAF.
For space reasons, TIFFs and Masks are not included.
The `processMNQuant.R` script will read these files and generate the plots in `Output/Plots/`

## Sankey Fate

An R-based script to generate Sankey diagrams from manual fate mapping.

## Notes

Mitotic progression figures were generated using [MitoticTiming](https://github.com/quantixed/MitoticTiming)