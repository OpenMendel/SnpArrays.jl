# BEDFiles.jl
Routines for reading and manipulating GWAS data in .bed files

## Background

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
The three possible types that can be observed are:
homozygous allele 1, coded as `0x00`, heterzygous, coded as `0x10`, and homozygous allele 2, coded as `0x11`.
Missing values are coded as `0x01`.

A single column - one SNP position over all `m` individuals - is packed into an
array of `div(m + 3, 4)` bytes (`UInt8` values).

## Installation

This package requires Julia v0.7.0 or later, which can be obtained from
<https://julialang.org/downloads/> or by building Julia from the sources in the
<https://github.com/JuliaLang/julia> repository.

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

## Loading a .bed file

The `BEDFile` struct contains the read-only, memory-mapped `.bed` file as a `Matrix{UInt8}`,
along with `m`, the number of individuals.
```@docs
BEDFile
```
For convenience, two `Int` matrices, `columncounts` and `rowcounts` are allocated but not populated until used.

The columns correspond to SNP positions.
Rows of the internal matrix are packed values from groups of 4 individuals.

```julia
julia> using BenchmarkTools, BEDFiles

julia> const bf = BEDFile(BEDFiles.datadir("mouse.bed"));

julia> size(bf)      # the virtual size of the GWAS data - 1940 observations at each of 10150 SNP positions
(1940, 10150)

julia> size(bf.data) # the actual size of the memory-mapped matrix of UInt8s
(485, 10150)
```

As described above, a column, consisting of `m` values in the range `0x00` to `0x03`, is packed into `div(m + 3, 4)` bytes.

(The calculation `div(i + 3, 4)` or `(i + 3) ÷ 4` occurs in some important loops in the code where it is evaluated as the
equivalent, but somewhat faster, expression `(i + 3) >> 2` that performs the integer division by 4 via shifting the number
two bits to the right.)

The virtual number of rows, `m`, can be given as a second argument in the call to `BEDFile`.
If omitted, `m` is determined as the number of lines in the `.fam` file. 

Because the file is memory-mapped this operation is fast, even for very large `.bed` files.
```julia
julia> @benchmark(BEDFile(BEDFiles.datadir("mouse.bed")))
BenchmarkTools.Trial: 
  memory estimate:  390.42 KiB
  allocs estimate:  82
  --------------
  minimum time:     127.349 μs (0.00% GC)
  median time:      135.760 μs (0.00% GC)
  mean time:        150.356 μs (7.05% GC)
  maximum time:     41.823 ms (99.33% GC)
  --------------
  samples:          10000
  evals/sample:     1
```

This file, from a study published in 2006, is about 5 Mb in size but data from recent studies, which have samples from tens of
thousands of individuals at over a million SNP positions, would be in the tens or even hundreds of Gb range.
## Raw summaries

Counts of each the four possible values for each column are returned by `counts`.

```julia
julia> counts(bf, dims=1)
4×10150 Array{Int64,2}:
  358   359  252   358    33   359    33  186   360  …    53    56    56    56    56    56    56    56
    2     0    4     3     4     1     4    1     3      171   174   173   173   162   173   174   175
 1003  1004  888  1004   442  1004   481  803  1002      186   242   242   242   242   242   242   242
  577   577  796   575  1461   576  1422  950   575     1530  1468  1469  1469  1480  1469  1468  1467
```

Column 2 has no missing values (code `0x01`, the second row in the column-counts table).
In that SNP position for this sample, 359 indivduals are homozygous allele 1 (`G` according to the `.bim` file), 1004 are heterozygous,
and 577 are homozygous allele 2 (`A`).

The counts by column and by row are cached in the `BEDFile` object.
Accesses after the first are extremely fast.
```julia
julia> @benchmark counts($bf, dims=1)
BenchmarkTools.Trial: 
  memory estimate:  0 bytes
  allocs estimate:  0
  --------------
  minimum time:     5.115 ns (0.00% GC)
  median time:      5.298 ns (0.00% GC)
  mean time:        5.244 ns (0.00% GC)
  maximum time:     23.231 ns (0.00% GC)
  --------------
  samples:          10000
  evals/sample:     1000
```

## Instantiating as a count of the second allele

In some operations on GWAS data the data are converted to counts of the second allele, according to

|BEDFile|count   |
|------:|--------:|
| 0x00  | 0       |
| 0x01  | missing |
| 0x10  | 1       |
| 0x11  | 2       |

This can be accomplished by indexing `bedvals`
```@docs
bedvals
```
with the `BEDFile` or with a view of the `BEDFile`,
producing an array of type `Union{Missing,Int8}`, which is the preferred way in v0.7 of
representing arrays that may contain missing values.
```julia
julia> bedvals[bf]
1940×10150 Array{Union{Missing, Int8},2}:
 1  1  1  1  2  1  2  1  1  1  1  2  1  1  2  …  2         2         2       
 1  1  2  1  1  1  1  2  1  1  1  1  2  2  1     2         2         2       
 2  2  2  2  2  2  2  2  2  2  2  2  2  2  2     2         2         2       
 1  1  1  1  1  1  1  2  1  1  1  1  2  2  1     2         2         2       
 2  2  2  2  2  2  2  2  2  2  2  2  2  2  2     1         1         1       
 1  1  1  1  2  1  2  1  1  1  1  2  1  1  2  …  2         2         2       
 1  1  1  1  2  1  2  1  1  1  1  2  1  1  2     2         2         2       
 1  1  2  1  1  1  1  2  1  1  1  1  2  2  1     2         2         2       
 1  1  2  1  1  1  1  2  1  1  1  1  2  2  1     2         2         2       
 ⋮              ⋮              ⋮              ⋱                              
 1  1  1  1  2  1  2  1  1  1  1  2  1  1  2     2         2         2       
 2  2  2  2  2  2  2  2  2  2  2  2  2  2  2     2         2         2       
 1  1  1  1  2  1  2  1  1  1  1  2  1  1  2     2         2         2       
 2  2  2  2  2  2  2  2  2  2  2  2  2  2  2     2         2         2       
 1  1  1  1  2  1  2  1  1  1  1  2  1  1  2  …  2         2         2       
 1  1  2  1  1  1  1  2  1  1  1  1  2  2  1     2         2         2       
 1  1  2  1  1  1  1  2  1  1  1  1  2  2  1     2         2         2       
 1  1  1  1  1  1  1  2  1  1  1  1  2  2  1      missing   missing   missing
 0  0  0  0  2  0  2  0  0  0  0  2  0  0  2     2         2         2       

julia> sort(unique(ans))
4-element Array{Union{Missing, Int8},1}:
 0       
 1       
 2       
  missing
```

### Summary statistics

The package provides methods for the generics `mean` and `var` from the `Statistics` package.
```julia
julia> mean(bf, dims=1)
1×10150 Array{Float64,2}:
 1.113  1.11237  1.28099  1.11203  …  1.8009  1.79966  1.79955  1.79943

julia> var(bf, dims=1)
1×10150 Array{Float64,2}:
 0.469929  0.470089  0.462605  0.469365  …  0.223714  0.223818  0.223923
```

These methods make use of the cached column or row counts and thus are very fast
```julia
julia> @benchmark mean(bf, dims=1)
BenchmarkTools.Trial: 
  memory estimate:  79.39 KiB
  allocs estimate:  2
  --------------
  minimum time:     37.777 μs (0.00% GC)
  median time:      38.186 μs (0.00% GC)
  mean time:        45.207 μs (14.53% GC)
  maximum time:     43.815 ms (99.83% GC)
  --------------
  samples:          10000
  evals/sample:     1
```

The column-wise or row-wise standard deviations are returned by `std`.
```julia
julia> std(bf, dims=2)
1940×1 Array{Float64,2}:
 0.6504997290784408
 0.6379008244533891
 0.6558172726141286
 0.6532675479248437
 0.6744432174014563
 0.6519092298111158
 0.6779881845456428
 0.6955814098050999
 0.6437566832989493
 ⋮                 
 0.6613451651895259
 0.6659810347614777
 0.6274577846909379
 0.6823658517777204
 0.6695299551061924
 0.710756592739754 
 0.6387913736114869
 0.6736492722732016
 0.688855476425891 
```

## Location of the missing values


The positions of the missing data are evaluated by
```@docs
missingpos
```
```julia
julia> mp = missingpos(bf)
1940×10150 SparseArrays.SparseMatrixCSC{Bool,Int32} with 33922 stored entries:
  [702  ,     1]  =  true
  [949  ,     1]  =  true
  [914  ,     3]  =  true
  [949  ,     3]  =  true
  [1604 ,     3]  =  true
  [1891 ,     3]  =  true
  [81   ,     4]  =  true
  [990  ,     4]  =  true
  [1882 ,     4]  =  true
  ⋮
  [1848 , 10150]  =  true
  [1851 , 10150]  =  true
  [1853 , 10150]  =  true
  [1860 , 10150]  =  true
  [1873 , 10150]  =  true
  [1886 , 10150]  =  true
  [1894 , 10150]  =  true
  [1897 , 10150]  =  true
  [1939 , 10150]  =  true
julia> @benchmark missingpos($bf)
BenchmarkTools.Trial: 
  memory estimate:  1.81 MiB
  allocs estimate:  19273
  --------------
  minimum time:     38.009 ms (0.00% GC)
  median time:      38.161 ms (0.00% GC)
  mean time:        38.761 ms (1.25% GC)
  maximum time:     80.250 ms (52.34% GC)
  --------------
  samples:          129
  evals/sample:     1
```

So, for example, the number of missing data values in each column can be evaluated as
```julia
julia> sum(mp, dims=1)
1×10150 Array{Int64,2}:
 2  0  4  3  4  1  4  1  3  3  0  4  0  …  174  173  173  162  173  174  175
```
although it is faster, but somewhat more obscure, to use
```
julia> view(counts(bf, dims=1), 2:2, :)
1×10150 view(::Array{Int64,2}, 2:2, :) with eltype Int64:
 2  0  4  3  4  1  4  1  3  3  0  4  0  …  174  173  173  162  173  174  175
```
