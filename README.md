# SnpArrays.jl

| **Documentation** | **Build Status** | **Code Coverage**  |
|-------------------|------------------|--------------------|
| [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://OpenMendel.github.io/SnpArrays.jl/latest) [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://OpenMendel.github.io/SnpArrays.jl/stable) | [![Build Status](https://travis-ci.org/OpenMendel/SnpArrays.jl.svg?branch=master)](https://travis-ci.org/OpenMendel/SnpArrays.jl) [![Build status](https://ci.appveyor.com/api/projects/status/wvaqu7i3ty2gk377/branch/master?svg=true)](https://ci.appveyor.com/project/Hua-Zhou/snparrays-jl-qavxa/branch/master) | [![Coverage Status](https://coveralls.io/repos/github/OpenMendel/SnpArrays.jl/badge.svg?branch=master)](https://coveralls.io/github/OpenMendel/SnpArrays.jl?branch=master) [![codecov](https://codecov.io/gh/OpenMendel/SnpArrays.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/OpenMendel/SnpArrays.jl) |

Routines for reading and manipulating compressed storage of biallelic [SNP (Single-Nucleotide Polymorphism)](https://en.wikipedia.org/wiki/Single-nucleotide_polymorphism) data.

Data from [*genome-wide association studies (GWAS)*](https://en.wikipedia.org/wiki/Genome-wide_association_study) are often saved as a [**PLINK binary biallelic genotype table**](https://www.cog-genomics.org/plink2/formats#bed) or `.bed` file.
To be useful, such files should be accompanied by a `.fam` file, containing metadata on the rows of the table, and a `.bim` file, containing metadata on the columns. The `.fam` and `.bim` files are in tab-separated format.

## Installation

This package requires Julia v0.7 or later, which can be obtained from
https://julialang.org/downloads/ or by building Julia from the sources in the
https://github.com/JuliaLang/julia repository.

The package has not yet been registered and must be installed using the repository location.
Start julia and use the `]` key to switch to the package manager REPL
```julia
(v1.1) pkg> add https://github.com/OpenMendel/SnpArrays.jl
```
Use the backspace key to return to the Julia REPL.

## Citation

If you use [OpenMendel](https://openmendel.github.io) analysis packages in your research, please cite the following reference in the resulting publications:

*Lange K, Papp JC, Sinsheimer JS, Sripracha R, Zhou H, Sobel EM (2013) Mendel: The Swiss army knife of genetic analysis programs. Bioinformatics 29:1568-1570.*

## Acknowledgments

Current implementation incorporates ideas in the package [BEDFiles.jl](https://github.com/dmbates/BEDFiles.jl) by Doug Bates (@dmbates).

This project is supported by the National Institutes of Health under NIGMS awards R01GM053275 and R25GM103774 and NHGRI award R01HG006139.
