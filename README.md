This code is actively developed at https://gitlab.umich.edu/VeatchLab/smlm-spacetime-resolution 

# Veatch Lab single molecule localization microscopy analysis software

This repository contains a variety of MATLAB code developed in the
[Veatch Lab](https://sites.lsa.umich.edu/veatch-lab/), which we use to
analyze point data from single molecule localization microscopy. These tools
are mostly based on pair correlation functions of various kinds, including
spatiotemporal auto- and cross-correlations. We apply these tools to
biophysical questions such as weak co-clustering of plasma membrane proteins,
diffusion of plasma membrane proteins, and characterization of the effective
precision, or "localization spread function" of the SMLM analysis pipeline.

We have developed these tools largely in tandem with two manuscripts:

TR Shaw, FJ Fazekas, S Kim, JC Flanagan-Natoli, ER Sumrall, SL Veatch. 2022.
Estimating the localization spread function of static single-molecule
localization microscopy images. *Biophys J* 121, 2906-2920. 
[doi:10.1016/j.bpj.2022.06.036](https://doi.org/10.1016/j.bpj.2022.06.036)

SA Shelby, TR Shaw, SL Veatch.
Measuring the co-localization and dynamics of mobile proteins in live cells
undergoing signaling responses. 2022. In preparation.

## Dependencies

The following MATLAB toolboxes are used by this software:
- Statistics and Machine Learning Toolbox
- Mapping Toolbox
- Curve Fitting Toolbox
- Image Processing Toolbox (only required by `spacewin_gui`)

## Installation

1. Add src/ and all of its subdirectories to your matlab path
```
addpath(genpath('src/'))
```
2. Compile the necessary mex files in `src/closepairs`
```
cd src/closepairs
compile_closepairs
```

## Brief description of contents

Pair correlation function calculations in several variants:
```
% correlation functions in space and time
spacetime_acor % autocorrelations (i.e. correlations between different points of the same type), as a function of distance and time separation
spacetime_acor_xy % autocorrelations as a function of x,y separation and time separation
spacetime_xcor % cross-correlations (i.e. between points of different types), as a function of distance and time separation

% correlation functions in space only
spatial_acor % autocorrelation in distance and time separation
spatial_xcor % cross-correlation in distance and time separation
```

Functions to compute effective "localization spread function"s (LSFs).
```
spacetime_resolution % LSF computation with angular average, i.e. as function of distance and time separation
spacetime_lsf_2d % LSF computation without angular average, i.e. as a function of x,y separation and time separation
```

Utility functions for dealing with point data:
```
spacewin_gui % A graphical user interface for drawing spatial windows (regions of interest)
default_iref % utility function to produce a good default imref2d on which to reconstruct a point dataset
reconstruct % generate a reconstructed image from a point dataset
gausblur % convolve an image with a 2d gaussian blur
```

Utility functions related to spatial and temporal windows (regions of interest)
are included under `src/windows/`, and utility functions for efficiently extracting
pairs of nearby localizations (in some combination of space and time) are
included in `src/crosspairs/`.

More detailed help text is provided with each function and can be referenced
with MATLAB's `help`.

## Examples

There is an example LSF estimation script in `src/example/example_NPC.m`, 
based on the NPC dataset from Shaw et al 2022.  Follow along with that to learn
how to use the LSF estimation code.

Scripts to replicate three of the figures of Shelby et al 2022 are also included,
under `src/example/xcor-protocol-figs/`.
