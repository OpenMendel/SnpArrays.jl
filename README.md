# SnpArrays.jl

| **Documentation** | **Build Status** | **Code Coverage**  | **Citation**  |
|-------------------|------------------|--------------------|--------------------|
| [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://OpenMendel.github.io/SnpArrays.jl/latest) [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://OpenMendel.github.io/SnpArrays.jl/stable) | [![build Actions Status](https://github.com/OpenMendel/SnpArrays.jl/workflows/CI/badge.svg)](https://github.com/OpenMendel/SnpArrays.jl/actions) | [![Coverage Status](https://coveralls.io/repos/github/OpenMendel/SnpArrays.jl/badge.svg?branch=master)](https://coveralls.io/github/OpenMendel/SnpArrays.jl?branch=master) [![codecov](https://codecov.io/gh/OpenMendel/SnpArrays.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/OpenMendel/SnpArrays.jl) | [![DOI](https://zenodo.org/badge/59141808.svg)](https://zenodo.org/badge/latestdoi/59141808) |

Routines for reading and manipulating compressed storage of biallelic [SNP (Single-Nucleotide Polymorphism)](https://en.wikipedia.org/wiki/Single-nucleotide_polymorphism) data.

Data from [*genome-wide association studies (GWAS)*](https://en.wikipedia.org/wiki/Genome-wide_association_study) are often saved as a [**PLINK binary biallelic genotype table**](https://www.cog-genomics.org/plink2/formats#bed) or `.bed` file.
To be useful, such files should be accompanied by a `.fam` file, containing metadata on the rows of the table, and a `.bim` file, containing metadata on the columns. The `.fam` and `.bim` files are in tab-separated format.

Linear algebra operations on the PLINK formatted data now support multi-threading and GPU (CUDA) computing.

## Installation

This package requires Julia v1.5 or later, which can be obtained from
https://julialang.org/downloads/ or by building Julia from the sources in the
https://github.com/JuliaLang/julia repository.

This package is registered in the default Julia package registry, and can be installed through standard package installation procedure.
```julia
using Pkg
pkg"add SnpArrays"
```
Use the backspace key to return to the Julia REPL.

## Citation

If you use [OpenMendel](https://openmendel.github.io) analysis packages in your research, please cite the following reference in the resulting publications:

*Zhou H, Sinsheimer JS, Bates DM, Chu BB, German CA, Ji SS, Keys KL, Kim J, Ko S, Mosher GD, Papp JC, Sobel EM, Zhai J, Zhou JJ, Lange K. OPENMENDEL: a cooperative programming project for statistical genetics. Hum Genet. 2020 Jan;139(1):61-71. doi: 10.1007/s00439-019-02001-z. Epub 2019 Mar 26. PMID: 30915546; PMCID: [PMC6763373](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6763373/).*

## Acknowledgments

Current implementation incorporates ideas in the package [BEDFiles.jl](https://github.com/dmbates/BEDFiles.jl) by Doug Bates (@dmbates).

Chris Elrod (@chriselrod) helped us accelerate CPU linear algebra through his great support of [LoopVectorization.jl](https://github.com/chriselrod/LoopVectorization.jl) package.

This project has been supported by the National Institutes of Health under awards R01GM053275, R01HG006139, R25GM103774, and 1R25HG011845.
