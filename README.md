# Installation

1. Add src/ and all of its subdirectories to your matlab path
```
addpath(genpath('src/'))
```
2. Compile the necessary mex files in `src/closepairs`
```
cd src/closepairs
compile_closepairs
```

# Method

The PSF estimation method is described and validated in a paper at
[arXiv:2202.08798 (physics.bio-ph)](https://arxiv.org/abs/2202.08798)
# Example

There is an example script in `src/example/example_NPC.m`, which applies the
PSF estimation method to the NPC dataset from the paper. Follow along with that to
learn how to use the resolution estimation code.

