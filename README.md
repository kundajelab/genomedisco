# genomedisco
[![Build Status](https://travis-ci.org/kundajelab/genomedisco.svg?branch=master)](https://travis-ci.org/kundajelab/genomedisco)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/kundajelab/genomedisco/blob/master/LICENSE)

`genomedisco` is a package for comparing contact maps of 3D genome structures, obtained from experiments such as Hi-C, Capture-C, ChIA-PET, HiChip, etc. It uses graph diffusion to smooth contact maps, and then compares them, resulting in a reproducibility score that can be used for quality control of biological replicates.

Installation
---

1. Install [Anaconda](https://www.continuum.io/downloads). 
2. Install genomedisco with the following command:
```
conda install genomedisco -c kundajelab
```
genomedisco is compatible with Python 2.

Quick start
---
To start, you will need the following:
- 2 contact maps for which you want to compute reproducibility


