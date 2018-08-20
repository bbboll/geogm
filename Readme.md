# GeoGM: Image labeling with graphical models
This is an implementation of the algorithm proposed in the paper [Image labeling based on graphical models using Wasserstein messages and geometric assignment](https://arxiv.org/abs/1710.01493).

## Prerequisites
Building is intended to run on Mac or Linux with clang. In both cases, the boost libraries (tested on version 1.67.01) are required.
On Linux, you also need to install [openBLAS](https://github.com/xianyi/OpenBLAS). The install location can be modified in the Makefile.