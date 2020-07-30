# SnpArrays.jl

| **Documentation** | **Build Status** | **Code Coverage**  | **Citation**  |
|-------------------|------------------|--------------------|--------------------|
| [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://OpenMendel.github.io/SnpArrays.jl/latest) [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://OpenMendel.github.io/SnpArrays.jl/stable) | [![Build Status](https://travis-ci.org/OpenMendel/SnpArrays.jl.svg?branch=master)](https://travis-ci.org/OpenMendel/SnpArrays.jl) [![Build status](https://ci.appveyor.com/api/projects/status/wvaqu7i3ty2gk377/branch/master?svg=true)](https://ci.appveyor.com/project/Hua-Zhou/snparrays-jl-qavxa/branch/master) | [![Coverage Status](https://coveralls.io/repos/github/OpenMendel/SnpArrays.jl/badge.svg?branch=master)](https://coveralls.io/github/OpenMendel/SnpArrays.jl?branch=master) [![codecov](https://codecov.io/gh/OpenMendel/SnpArrays.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/OpenMendel/SnpArrays.jl) | [![DOI](https://zenodo.org/badge/59141808.svg)](https://zenodo.org/badge/latestdoi/59141808) |

Routines for reading and manipulating compressed storage of biallelic [SNP (Single-Nucleotide Polymorphism)](https://en.wikipedia.org/wiki/Single-nucleotide_polymorphism) data.

Data from [*genome-wide association studies (GWAS)*](https://en.wikipedia.org/wiki/Genome-wide_association_study) are often saved as a [**PLINK binary biallelic genotype table**](https://www.cog-genomics.org/plink2/formats#bed) or `.bed` file.
To be useful, such files should be accompanied by a `.fam` file, containing metadata on the rows of the table, and a `.bim` file, containing metadata on the columns. The `.fam` and `.bim` files are in tab-separated format.

## Installation

This package requires Julia v1.4 or later, which can be obtained from
https://julialang.org/downloads/ or by building Julia from the sources in the
https://github.com/JuliaLang/julia repository.

The package has not yet been registered and must be installed using the repository location.
Start julia and use the `]` key to switch to the package manager REPL
```julia
(v1.4) pkg> add https://github.com/OpenMendel/SnpArrays.jl
```
Use the backspace key to return to the Julia REPL.

## Citation

If you use [OpenMendel](https://openmendel.github.io) analysis packages in your research, please cite the following reference in the resulting publications:

*OPENMENDEL: a cooperative programming project for statistical genetics. Zhou H, Sinsheimer JS, Bates DM, Chu BB, German CA, Ji SS, Keys KL, Kim J, Ko S, Mosher GD, Papp JC, Sobel EM, Zhai J, Zhou JJ, Lange K. Hum Genet. 2019 Mar 26. doi: 10.1007/s00439-019-02001-z. [Epub ahead of print] PMID: 30915546*

## Acknowledgments

Current implementation incorporates ideas in the package [BEDFiles.jl](https://github.com/dmbates/BEDFiles.jl) by Doug Bates (@dmbates).

GitHub user @chriselrod helped us accelerate CPU linear algebra through his great support of [LoopVectorization.jl](https://github.com/chriselrod/LoopVectorization.jl) package.

This project is supported by the National Institutes of Health under NIGMS awards R01GM053275 and R25GM103774 and NHGRI award R01HG006139.
