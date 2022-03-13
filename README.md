# SMLM spacetime resolution

This software provides the computations described in

TR Shaw, FJ Fazekas, S Kim, JC Flanagan-Natoli, ER Sumrall, SL Veatch.
A method to estimate the effective point spread function of static single
molecule localization microscopy images. Biorxiv. 2022.03.05.483117;
[doi:10.1101/2022.03.05.483117](https://www.biorxiv.org/content/10.1101/2022.03.05.483117v1).

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

## Example

There is an example script in `src/example/example_NPC.m`, which applies the
PSF estimation method to the NPC dataset from the paper. Follow along with that to
learn how to use the resolution estimation code.

