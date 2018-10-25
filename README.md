# BEDFiles.jl
Routines for reading and manipulating GWAS data in .bed files

| **Documentation**                                                               | **PackageEvaluator**                                            | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:---------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-latest-img]][docs-latest-url] | [![][pkg-0.7-img]][pkg-0.7-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] [![][coveralls-img]][coveralls-url] [![][codecov-img]][codecov-url]|

Data from [*Genome-wide association studies*](https://en.wikipedia.org/wiki/Genome-wide_association_study)
are often saved as a [**PLINK binary biallelic genotype table**](https://www.cog-genomics.org/plink2/formats#bed)
or `.bed` file.
To be useful, such files should be accompanied by a `.fam` file, containing metadata on the rows of the table, and a `.bim` file,
containing metadata on the columns.
The `.fam` and `.bim` files are in tab-separated format.

The table contains the observed allelic type at `n`
[*single-nucleotide polymorphism*](https://en.wikipedia.org/wiki/Single-nucleotide_polymorphism) (SNP) positions 
for `m` individuals.

A SNP corresponds to a nucleotide position on the genome where some degree of variation has been observed in a population,
with each individual have one of two possible *alleles* at that position on each of a pair of chromosomes.
Three possible types can be observed are:
homozygous allele 1, coded as `0x00`, heterzygous, coded as `0x10`, and homozygous allele 2, coded as `0x11`.
Missing values are coded as `0x01`.

A single column - one SNP position over all `m` individuals - is packed into an
array of `div(m + 3, 4)` bytes (`UInt8` values).

## Installation

This package requires Julia v0.7.0 or later, which can be obtained from
https://julialang.org/downloads/ or by building Julia from the sources in the
https://github.com/JuliaLang/julia repository.

The package has not yet been registered and must be installed using the repository location.
Start julia and use the `]` key to switch to the package manager REPL
```julia
(v0.7) pkg> add https://github.com/dmbates/BEDFiles.jl.git#master
  Updating git-repo `https://github.com/dmbates/BEDFiles.jl.git`
  Updating registry at `~/.julia/registries/Uncurated`
  Updating git-repo `https://github.com/JuliaRegistries/Uncurated.git`
 Resolving package versions...
  Updating `~/.julia/environments/v0.7/Project.toml`
  [6f44c9a6] + BEDFiles v0.1.0 #master (https://github.com/dmbates/BEDFiles.jl.git)
  Updating `~/.julia/environments/v0.7/Manifest.toml`
  [6f44c9a6] + BEDFiles v0.1.0 #master (https://github.com/dmbates/BEDFiles.jl.git)
  [6fe1bfb0] + OffsetArrays v0.6.0
  [10745b16] + Statistics 
```

Use the backspace key to return to the Julia REPL.

Please see the documentation [![][docs-latest-img]][docs-latest-url] for usage.

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://dmbates.github.io/BEDFiles.jl/latest

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://dmbates.github.io/BEDFiles.jl/stable

[travis-img]: https://travis-ci.org/dmbates/BEDFiles.jl.svg?branch=master
[travis-url]: https://travis-ci.org/dmbates/BEDFiles.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/bifqhte27nekp97m?svg=true
[appveyor-url]: https://ci.appveyor.com/project/dmbates/bedfiles-jl/branch/master

[coveralls-img]: https://coveralls.io/repos/github/dmbates/BEDFiles.jl/badge.svg?branch=master
[coveralls-url]: https://coveralls.io/github/dmbates/BEDFiles.jl?branch=master

[codecov-img]: https://codecov.io/github/dmbates/BEDFiles.jl/badge.svg?branch=master
[codecov-url]: https://codecov.io/github/dmbates/BEDFiles.jl?branch=master

[issues-url]: https://github.com/dmbates/BEDFiles.jl/issues

[pkg-0.7-img]: http://pkg.julialang.org/badges/BEDFiles_0.7.svg
[pkg-0.7-url]: http://pkg.julialang.org/?pkg=BEDFiles
