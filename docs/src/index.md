
# SnpArrays.jl

Data from [*genome-wide association studies (GWAS)*](https://en.wikipedia.org/wiki/Genome-wide_association_study) are often saved as a [**PLINK binary biallelic genotype table**](https://www.cog-genomics.org/plink2/formats#bed) or `.bed` file. To be useful, such files should be accompanied by a `.fam` file, containing metadata on the rows of the table, and a `.bim` file,
containing metadata on the columns. The `.fam` and `.bim` files are in tab-separated format.

The table contains the observed allelic type at `n` [*single nucleotide polymorphism*](https://en.wikipedia.org/wiki/Single-nucleotide_polymorphism) (SNP) positions for `m` individuals. A SNP corresponds to a nucleotide position on the genome where some degree of variation has been observed in a population, with each individual have one of two possible *alleles* at that position on each of a pair of chromosomes. Three possible genotypes and corresponding coding are

| Genotype | Plink/SnpArray |  
|:---:|:---:|  
| A1,A1 | 0x00 |  
| missing | 0x01 |
| A1,A2 | 0x02 |  
| A2,A2 | 0x03 |  

## Installation

This package requires Julia v0.7 or later, which can be obtained from
<https://julialang.org/downloads/> or by building Julia from the sources in the
<https://github.com/JuliaLang/julia> repository.

The package has not yet been registered and must be installed using the repository location.
Start julia and use the `]` key to switch to the package manager REPL
```julia
(v1.1) pkg> add https://github.com/OpenMendel/SnpArrays.jl
```
Use the backspace key to return to the Julia REPL.


```julia
versioninfo()
```

    Julia Version 1.1.0
    Commit 80516ca202 (2019-01-21 21:24 UTC)
    Platform Info:
      OS: macOS (x86_64-apple-darwin14.5.0)
      CPU: Intel(R) Core(TM) i7-6920HQ CPU @ 2.90GHz
      WORD_SIZE: 64
      LIBM: libopenlibm
      LLVM: libLLVM-6.0.1 (ORCJIT, skylake)
    Environment:
      JULIA_EDITOR = code



```julia
# for use in this tutorial
using SnpArrays, BenchmarkTools, Glob
```

## Example data

There are two example data sets attached to this package. They are availabe in the `data` folder of the package.


```julia
datapath = normpath(SnpArrays.datadir())
```




    "/Users/huazhou/.julia/dev/SnpArrays/data"




```julia
readdir(glob"mouse.*", datapath)
```




    3-element Array{String,1}:
     "/Users/huazhou/.julia/dev/SnpArrays/data/mouse.bed"
     "/Users/huazhou/.julia/dev/SnpArrays/data/mouse.bim"
     "/Users/huazhou/.julia/dev/SnpArrays/data/mouse.fam"



Data set `EUR_subset` contains no missing genotypes. It is located at


```julia
readdir(glob"EUR_subset.*", datapath)
```




    3-element Array{String,1}:
     "/Users/huazhou/.julia/dev/SnpArrays/data/EUR_subset.bed"
     "/Users/huazhou/.julia/dev/SnpArrays/data/EUR_subset.bim"
     "/Users/huazhou/.julia/dev/SnpArrays/data/EUR_subset.fam"



Data from recent studies, which have samples from tens of thousands of individuals at over a million SNP positions, would be in the tens or even hundreds of Gb range.

# SnpArray

`SnpArray` is the fundamental type for dealing with genotype data in Plink bed file. Each row of `SnpArray` is a sample and each column a SNP.

## Constructor

There are various ways to initialize a SnpArray.

### Intitialize from Plink file set

SnpArray can be initialized from the Plink bed file. The corresponding `.fam` needs to be present, which is used to determine the number of individuals.


```julia
const mouse = SnpArray(SnpArrays.datadir("mouse.bed"))
```




    1940×10150 SnpArray:
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x02  0x02  0x02  0x02  0x02  0x02
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x02  0x02  0x02  0x02  0x02  0x02
     0x03  0x03  0x03  0x03  0x03  0x03  …  0x00  0x00  0x00  0x00  0x00  0x00
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x00  0x00  0x00  0x00  0x00  0x00
        ⋮                             ⋮  ⋱           ⋮                        
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02     0x01  0x01  0x01  0x01  0x01  0x01
     0x00  0x00  0x00  0x00  0x03  0x00     0x03  0x03  0x03  0x03  0x03  0x03



The virtual size of the GWAS data is 1940 observations at each of 10150 SNP positions.


```julia
size(mouse)
```




    (1940, 10150)



Because the file is memory-mapped opening the file and accessing the data is fast, even for very large .bed files.


```julia
@btime(SnpArray(SnpArrays.datadir("mouse.bed")));
```

      115.086 μs (90 allocations: 391.09 KiB)


By default, the memory-mapped file is read only, changing entries is not allowed.


```julia
mouse[1, 1] = 0x00
```


    ReadOnlyMemoryError()

    

    Stacktrace:

     [1] | at ./int.jl:320 [inlined]

     [2] setindex!(::SnpArray, ::UInt8, ::Int64, ::Int64) at /Users/huazhou/.julia/dev/SnpArrays/src/snparray.jl:131

     [3] top-level scope at In[9]:1


To possibly change genoytpes in a bed file, open with write permission
```julia
mouse = SnpArray(SnpArrays.datadir("mouse.bed"), "w")
```

### Initialize from only bed file

If only the bed file is present, user is required to supply the number of individuals in the second argument.


```julia
SnpArray(SnpArrays.datadir("mouse.bed"), 1940)
```




    1940×10150 SnpArray:
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x02  0x02  0x02  0x02  0x02  0x02
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x02  0x02  0x02  0x02  0x02  0x02
     0x03  0x03  0x03  0x03  0x03  0x03  …  0x00  0x00  0x00  0x00  0x00  0x00
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x00  0x00  0x00  0x00  0x00  0x00
        ⋮                             ⋮  ⋱           ⋮                        
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02     0x01  0x01  0x01  0x01  0x01  0x01
     0x00  0x00  0x00  0x00  0x03  0x00     0x03  0x03  0x03  0x03  0x03  0x03



### Initialize from compressed Plink files

SnpArray can be initialized from Plink files in compressed formats: `gz`, `zlib`, `zz`, `xz`, `zst`, or `bz2`. For a complete list type
```julia
SnpArrays.ALLOWED_FORMAT
```
If you want to support a new compressed format, file an issue.

Let us first compress the mouse data in gz format. We see gz format takes less than 1/3 storage of original Plink files.


```julia
compress_plink(SnpArrays.datadir("mouse"), "gz")
readdir(glob"mouse.*.gz", datapath)
```




    3-element Array{String,1}:
     "/Users/huazhou/.julia/dev/SnpArrays/data/mouse.bed.gz"
     "/Users/huazhou/.julia/dev/SnpArrays/data/mouse.bim.gz"
     "/Users/huazhou/.julia/dev/SnpArrays/data/mouse.fam.gz"



To initialize SnpArray from gzipped Plink file, simply used the bed file with name ending with `.bed.gz`:


```julia
# requires corresponding `.fam.gz` file
SnpArray(SnpArrays.datadir("mouse.bed.gz"))
```




    1940×10150 SnpArray:
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x02  0x02  0x02  0x02  0x02  0x02
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x02  0x02  0x02  0x02  0x02  0x02
     0x03  0x03  0x03  0x03  0x03  0x03  …  0x00  0x00  0x00  0x00  0x00  0x00
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x00  0x00  0x00  0x00  0x00  0x00
        ⋮                             ⋮  ⋱           ⋮                        
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02     0x01  0x01  0x01  0x01  0x01  0x01
     0x00  0x00  0x00  0x00  0x03  0x00     0x03  0x03  0x03  0x03  0x03  0x03



or


```julia
# does not require corresponding `.fam.gz` file
SnpArray(SnpArrays.datadir("mouse.bed.gz"), 1940)
```




    1940×10150 SnpArray:
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x02  0x02  0x02  0x02  0x02  0x02
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x02  0x02  0x02  0x02  0x02  0x02
     0x03  0x03  0x03  0x03  0x03  0x03  …  0x00  0x00  0x00  0x00  0x00  0x00
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x00  0x00  0x00  0x00  0x00  0x00
        ⋮                             ⋮  ⋱           ⋮                        
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02     0x01  0x01  0x01  0x01  0x01  0x01
     0x00  0x00  0x00  0x00  0x03  0x00     0x03  0x03  0x03  0x03  0x03  0x03




```julia
# clean up
rm(SnpArrays.datadir("mouse.bed.gz"), force=true)
rm(SnpArrays.datadir("mouse.fam.gz"), force=true)
rm(SnpArrays.datadir("mouse.bim.gz"), force=true)
```

### Initialize and create bed file

Initialize 5 rows and 3 columns with all (A1, A1) genotype (0x00) and memory-map to a bed file `tmp.bed`


```julia
tmpbf = SnpArray("tmp.bed", 5, 3)
```




    5×3 SnpArray:
     0x00  0x00  0x00
     0x00  0x00  0x00
     0x00  0x00  0x00
     0x00  0x00  0x00
     0x00  0x00  0x00



Change entries


```julia
tmpbf[1:2, 1:2] .= 0x03
tmpbf
```




    5×3 SnpArray:
     0x03  0x03  0x00
     0x03  0x03  0x00
     0x00  0x00  0x00
     0x00  0x00  0x00
     0x00  0x00  0x00




```julia
fill!(tmpbf, 0x02)
tmpbf
```




    5×3 SnpArray:
     0x02  0x02  0x02
     0x02  0x02  0x02
     0x02  0x02  0x02
     0x02  0x02  0x02
     0x02  0x02  0x02




```julia
# clean up
rm("tmp.bed", force=true)
```

Initialize 5 rows and 3 columns with undefined genotypes without memory-mapping to any file


```julia
tmpbf = SnpArray(undef, 5, 3)
```




    5×3 SnpArray:
     0x00  0x00  0x00
     0x00  0x00  0x00
     0x00  0x00  0x00
     0x00  0x00  0x00
     0x00  0x00  0x00



Create a bed file corresponding to an existing SnpArray and memory-map it.


```julia
tmpbf = SnpArray("tmp.bed", tmpbf)
```




    5×3 SnpArray:
     0x00  0x00  0x00
     0x00  0x00  0x00
     0x00  0x00  0x00
     0x00  0x00  0x00
     0x00  0x00  0x00




```julia
tmpbf[1, 1] = 0x02
tmpbf
```




    5×3 SnpArray:
     0x02  0x00  0x00
     0x00  0x00  0x00
     0x00  0x00  0x00
     0x00  0x00  0x00
     0x00  0x00  0x00




```julia
# clean up
rm("tmp.bed", force=true)
```

## `convert` and `copyto!`

Most common usage of SnpArray is to convert genotypes to numeric values for statistical analysis. Conversion rule depends on genetic models (additive, dominant, or recessive), centering, scaling, or imputation.

### `convert`

`convert` function has 4 keyword arguments: `model`, `center`, `scale`, and `impute`.

`model` keyword specifies the SNP model for conversion. By default `convert` function translates genotypes according to the *additive* SNP model, which essentially counts the number of **A2** allele (0, 1 or 2) per genotype. Other SNP models are *dominant* and *recessive*, both in terms of the **A2** allele.

| Genotype | `SnpArray` | `model=ADDITIVE_MODEL` | `model=DOMINANT_MODEL` | `model=RECESSIVE_MODEL` |    
|:---:|:---:|:---:|:---:|:---:|  
| A1,A1 | 0x00 | 0 | 0 | 0 |  
| missing | 0x01 | NaN | NaN | NaN |
| A1,A2 | 0x02 | 1 | 1 | 0 |  
| A2,A2 | 0x03 | 2 | 1 | 1 |  

`center=true` tells `convert` to center each column by its mean. Default is `false`.

`scale=true` tells `convert` to scale each column by its standard deviation. Default is `false`.

`impute=true` tells `convert` to impute missing genotypes (0x01) by column mean. Default is `false`.

Convert whole SnpArray to a Float64 matrix using defaults (`model=ADDITIVE_MODEL`, `center=false`, `scale=false`, `impute=false`)


```julia
convert(Matrix{Float64}, mouse)
```




    1940×10150 Array{Float64,2}:
     1.0  1.0  1.0  1.0  2.0  1.0  2.0  1.0  …    2.0    2.0    2.0    2.0    2.0
     1.0  1.0  2.0  1.0  1.0  1.0  1.0  2.0       2.0    2.0    2.0    2.0    2.0
     2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0       2.0    2.0    2.0    2.0    2.0
     1.0  1.0  1.0  1.0  1.0  1.0  1.0  2.0       2.0    2.0    2.0    2.0    2.0
     2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0       1.0    1.0    1.0    1.0    1.0
     1.0  1.0  1.0  1.0  2.0  1.0  2.0  1.0  …    2.0    2.0    2.0    2.0    2.0
     1.0  1.0  1.0  1.0  2.0  1.0  2.0  1.0       2.0    2.0    2.0    2.0    2.0
     1.0  1.0  2.0  1.0  1.0  1.0  1.0  2.0       2.0    2.0    2.0    2.0    2.0
     1.0  1.0  2.0  1.0  1.0  1.0  1.0  2.0       2.0    2.0    2.0    2.0    2.0
     2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0       1.0    1.0    1.0    1.0    1.0
     2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0  …    0.0    0.0    0.0    0.0    0.0
     1.0  1.0  1.0  1.0  2.0  1.0  2.0  1.0       2.0    2.0    2.0    2.0    2.0
     2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0       0.0    0.0    0.0    0.0    0.0
     ⋮                        ⋮              ⋱    ⋮                              
     2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0       2.0    2.0    2.0    2.0    2.0
     2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0       2.0    2.0    2.0    2.0    2.0
     1.0  1.0  1.0  1.0  1.0  1.0  1.0  2.0  …    2.0    2.0    2.0    2.0    2.0
     1.0  1.0  1.0  1.0  2.0  1.0  2.0  1.0       2.0    2.0    2.0    2.0    2.0
     2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0       2.0    2.0    2.0    2.0    2.0
     1.0  1.0  1.0  1.0  2.0  1.0  2.0  1.0       2.0    2.0    2.0    2.0    2.0
     2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0       2.0    2.0    2.0    2.0    2.0
     1.0  1.0  1.0  1.0  2.0  1.0  2.0  1.0  …    2.0    2.0    2.0    2.0    2.0
     1.0  1.0  2.0  1.0  1.0  1.0  1.0  2.0       2.0    2.0    2.0    2.0    2.0
     1.0  1.0  2.0  1.0  1.0  1.0  1.0  2.0       2.0    2.0    2.0    2.0    2.0
     1.0  1.0  1.0  1.0  1.0  1.0  1.0  2.0     NaN    NaN    NaN    NaN    NaN  
     0.0  0.0  0.0  0.0  2.0  0.0  2.0  0.0       2.0    2.0    2.0    2.0    2.0



!!! note  

    When `convert` or `copyto!` a slice or subarray of SnpArray, using `view`, `@view` or `views` is necessary for both correctness and efficiency. Without view, it's simply converting the UInt8 coding in original bed file.
    

Convert a column to Float64 vector using defaults (`model=ADDITIVE_MODEL`, `center=false`, `scale=false`, `impute=false`).


```julia
# convert(Vector{Float64}, view(mouse, :, 1)) # alternative syntax
# @views convert(Vector{Float64}, mouse[:, 1]) # alternative syntax
convert(Vector{Float64}, @view(mouse[:, 1]))
```




    1940-element Array{Float64,1}:
     1.0
     1.0
     2.0
     1.0
     2.0
     1.0
     1.0
     1.0
     1.0
     2.0
     2.0
     1.0
     2.0
     ⋮  
     2.0
     2.0
     1.0
     1.0
     2.0
     1.0
     2.0
     1.0
     1.0
     1.0
     1.0
     0.0



Convert a subarray of SnpArray to Float64 matrix using defaults (`model=ADDITIVE_MODEL`, `center=false`, `scale=false`, `impute=false`).


```julia
convert(Matrix{Float64}, @view(mouse[1:2:10, 1:2:10]))
```




    5×5 Array{Float64,2}:
     1.0  1.0  2.0  2.0  1.0
     2.0  2.0  2.0  2.0  2.0
     2.0  2.0  2.0  2.0  2.0
     1.0  1.0  2.0  2.0  1.0
     1.0  2.0  1.0  1.0  1.0



Different SNP models (`ADDITIVE_MODEL` vs `DOMINANT_MODEL` vs `RECESSIVE_MODEL`)


```julia
@views [convert(Vector{Float64}, mouse[:, 1], model=ADDITIVE_MODEL) convert(Vector{Float64}, mouse[:, 1], model=DOMINANT_MODEL) convert(Vector{Float64}, mouse[:, 1], model=RECESSIVE_MODEL)]
```




    1940×3 Array{Float64,2}:
     1.0  1.0  0.0
     1.0  1.0  0.0
     2.0  1.0  1.0
     1.0  1.0  0.0
     2.0  1.0  1.0
     1.0  1.0  0.0
     1.0  1.0  0.0
     1.0  1.0  0.0
     1.0  1.0  0.0
     2.0  1.0  1.0
     2.0  1.0  1.0
     1.0  1.0  0.0
     2.0  1.0  1.0
     ⋮            
     2.0  1.0  1.0
     2.0  1.0  1.0
     1.0  1.0  0.0
     1.0  1.0  0.0
     2.0  1.0  1.0
     1.0  1.0  0.0
     2.0  1.0  1.0
     1.0  1.0  0.0
     1.0  1.0  0.0
     1.0  1.0  0.0
     1.0  1.0  0.0
     0.0  0.0  0.0



Center and scale (last column) while `convert`


```julia
convert(Vector{Float64}, @view(mouse[:, end]), center=true, scale=true)
```




    1940-element Array{Float64,1}:
       0.4721460731147541
       0.4721460731147541
       0.4721460731147541
       0.4721460731147541
      -1.8819155626127624
       0.4721460731147541
       0.4721460731147541
       0.4721460731147541
       0.4721460731147541
      -1.8819155626127624
      -4.2359771983402785
       0.4721460731147541
      -4.2359771983402785
       ⋮                 
       0.4721460731147541
       0.4721460731147541
       0.4721460731147541
       0.4721460731147541
       0.4721460731147541
       0.4721460731147541
       0.4721460731147541
       0.4721460731147541
       0.4721460731147541
       0.4721460731147541
     NaN                 
       0.4721460731147541



Center, scale, and impute (last column) while `convert`


```julia
convert(Vector{Float64}, @view(mouse[:, end]), center=true, scale=true, impute=true)
```




    1940-element Array{Float64,1}:
      0.4721460731147541
      0.4721460731147541
      0.4721460731147541
      0.4721460731147541
     -1.8819155626127624
      0.4721460731147541
      0.4721460731147541
      0.4721460731147541
      0.4721460731147541
     -1.8819155626127624
     -4.2359771983402785
      0.4721460731147541
     -4.2359771983402785
      ⋮                 
      0.4721460731147541
      0.4721460731147541
      0.4721460731147541
      0.4721460731147541
      0.4721460731147541
      0.4721460731147541
      0.4721460731147541
      0.4721460731147541
      0.4721460731147541
      0.4721460731147541
      0.0               
      0.4721460731147541



### `copyto!`

`copyto!` is the in-place version of `convert`. It takes the same keyword arguments (`model`, `center`, `scale`, `impute`) as `convert`.

Copy a column to a Float64 vector using defaults (`model=:additive`, `center=false`, `scale=false`, `impute=false`).


```julia
v = zeros(size(mouse, 1))
copyto!(v, @view(mouse[:, 1]))
```




    1940-element Array{Float64,1}:
     1.0
     1.0
     2.0
     1.0
     2.0
     1.0
     1.0
     1.0
     1.0
     2.0
     2.0
     1.0
     2.0
     ⋮  
     2.0
     2.0
     1.0
     1.0
     2.0
     1.0
     2.0
     1.0
     1.0
     1.0
     1.0
     0.0




```julia
@btime(copyto!($v, $@view(mouse[:, 1])));
```

      3.402 μs (0 allocations: 0 bytes)


Copy columns using defaults


```julia
v2 = zeros(size(mouse, 1), 2)
copyto!(v2, @view(mouse[:, 1:2]))
```




    1940×2 Array{Float64,2}:
     1.0  1.0
     1.0  1.0
     2.0  2.0
     1.0  1.0
     2.0  2.0
     1.0  1.0
     1.0  1.0
     1.0  1.0
     1.0  1.0
     2.0  2.0
     2.0  2.0
     1.0  1.0
     2.0  2.0
     ⋮       
     2.0  2.0
     2.0  2.0
     1.0  1.0
     1.0  1.0
     2.0  2.0
     1.0  1.0
     2.0  2.0
     1.0  1.0
     1.0  1.0
     1.0  1.0
     1.0  1.0
     0.0  0.0




```julia
# roughly double the cost of copying 1 column
@btime(copyto!($v2, $@view(mouse[:, 1:2])));
```

      6.879 μs (0 allocations: 0 bytes)


Center and scale


```julia
copyto!(v, @view(mouse[:, 1]), center=true, scale=true)
```




    1940-element Array{Float64,1}:
     -0.16084075452851265
     -0.16084075452851265
      1.2624897581484626 
     -0.16084075452851265
      1.2624897581484626 
     -0.16084075452851265
     -0.16084075452851265
     -0.16084075452851265
     -0.16084075452851265
      1.2624897581484626 
      1.2624897581484626 
     -0.16084075452851265
      1.2624897581484626 
      ⋮                  
      1.2624897581484626 
      1.2624897581484626 
     -0.16084075452851265
     -0.16084075452851265
      1.2624897581484626 
     -0.16084075452851265
      1.2624897581484626 
     -0.16084075452851265
     -0.16084075452851265
     -0.16084075452851265
     -0.16084075452851265
     -1.584171267205488  




```julia
# more cost becoz of extra pass for center, scale, and/or impute
@btime(copyto!($v, $(@view(mouse[:, 1])), center=true, scale=true));
```

      6.302 μs (0 allocations: 0 bytes)


Looping over all columns


```julia
v = Vector{Float64}(undef, size(mouse, 1))
function loop_test(v, s)
    for j in 1:size(s, 2)
        copyto!(v, @view(s[:, j]))
    end
end
@btime(loop_test($v, $mouse))
```

      41.091 ms (10150 allocations: 475.78 KiB)


Copy whole SnpArray


```julia
M = similar(mouse, Float64)
@btime(copyto!($M, $mouse));
```

      41.126 ms (0 allocations: 0 bytes)


## Summaries

### Counts

Counts of each the four possible values for each column are returned by `counts`.`


```julia
counts(mouse, dims=1)
```




    4×10150 Array{Int64,2}:
      358   359  252   358    33   359  …    56    56    56    56    56    56
        2     0    4     3     4     1      173   173   162   173   174   175
     1003  1004  888  1004   442  1004      242   242   242   242   242   242
      577   577  796   575  1461   576     1469  1469  1480  1469  1468  1467



Column 2 has no missing values (code `0x01`, the second row in the column-counts table).
In that SNP position for this sample, 359 indivduals are homozygous allele 1 (`G` according to the `.bim` file), 1004 are heterozygous, and 577 are homozygous allele 2 (`A`).

The counts by column and by row are cached in the `SnpArray` object. Accesses after the first are extremely fast.


```julia
@btime(counts($mouse, dims=1));
```

      6.792 ns (0 allocations: 0 bytes)


### Minor allele frequencies

Minor allele frequencies (MAF) for each SNP.


```julia
maf(mouse)
```




    10150-element Array{Float64,1}:
     0.4434984520123839  
     0.4438144329896907  
     0.359504132231405   
     0.4439855446566856  
     0.13119834710743805 
     0.44404332129963897 
     0.1412706611570248  
     0.30299123259412064 
     0.4445018069179143  
     0.44424367578729995 
     0.43427835051546393 
     0.14075413223140498 
     0.304639175257732   
     ⋮                   
     0.0527624309392265  
     0.052980132450331174
     0.08079096045197742 
     0.08253250423968339 
     0.08253250423968339 
     0.10022650056625138 
     0.10016977928692694 
     0.10016977928692694 
     0.09955005624296964 
     0.10016977928692694 
     0.10022650056625138 
     0.10028328611898019 



Minor allele (`false` means A1 is the minor allele; `true` means A2 is the minor allele) for each SNP.


```julia
minorallele(mouse)
```




    10150-element BitArray{1}:
     false
     false
     false
     false
     false
     false
     false
     false
     false
     false
     false
     false
     false
         ⋮
     false
     false
     false
     false
     false
     false
     false
     false
     false
     false
     false
     false



### `mean` and `var`

The package provides methods for the generics `mean` and `var` from the `Statistics` package.


```julia
mean(mouse, dims=1)
```




    1×10150 Array{Float64,2}:
     1.113  1.11237  1.28099  1.11203  …  1.8009  1.79966  1.79955  1.79943




```julia
mean(mouse, dims=1, model=DOMINANT_MODEL)
```




    1×10150 Array{Float64,2}:
     0.815273  0.814948  0.869835  0.815178  …  0.968308  0.96829  0.968272




```julia
var(mouse, dims=1)
```




    1×10150 Array{Float64,2}:
     0.469929  0.470089  0.462605  0.469365  …  0.223714  0.223818  0.223923



These methods make use of the cached column or row counts and thus are very fast


```julia
@btime(mean($mouse, dims=1));
```

      13.113 μs (2 allocations: 79.39 KiB)


The column-wise or row-wise standard deviations are returned by `std`.


```julia
std(mouse, dims=2)
```




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
     0.6505283141088536
     0.665444994623426 
     0.659392039592328 
     0.6641674726999468
     ⋮                 
     0.6599158250006595
     0.688387450736178 
     0.6664063015924304
     0.6613451651895259
     0.6659810347614777
     0.6274577846909379
     0.6823658517777204
     0.6695299551061924
     0.710756592739754 
     0.6387913736114869
     0.6736492722732016
     0.688855476425891 



### Missing rate

Proportion of missing genotypes


```julia
missingrate(mouse, 1)
```




    10150-element Array{Float64,1}:
     0.0010309278350515464
     0.0                  
     0.002061855670103093 
     0.0015463917525773195
     0.002061855670103093 
     0.0005154639175257732
     0.002061855670103093 
     0.0005154639175257732
     0.0015463917525773195
     0.0015463917525773195
     0.0                  
     0.002061855670103093 
     0.0                  
     ⋮                    
     0.06701030927835051  
     0.06597938144329897  
     0.08762886597938144  
     0.08814432989690722  
     0.08814432989690722  
     0.08969072164948454  
     0.08917525773195877  
     0.08917525773195877  
     0.08350515463917525  
     0.08917525773195877  
     0.08969072164948454  
     0.09020618556701031  




```julia
missingrate(mouse, 2)
```




    1940-element Array{Float64,1}:
     0.00019704433497536947
     0.0                   
     0.018423645320197045  
     0.0007881773399014779 
     0.0                   
     0.004236453201970443  
     0.0051231527093596055 
     0.00039408866995073894
     0.005517241379310344  
     0.0016748768472906405 
     0.0                   
     9.852216748768474e-5  
     0.0004926108374384236 
     ⋮                     
     0.000689655172413793  
     0.004729064039408867  
     0.0004926108374384236 
     0.001083743842364532  
     0.00019704433497536947
     0.0025615763546798028 
     0.0038423645320197044 
     0.001379310344827586  
     0.0064039408866995075 
     0.002857142857142857  
     0.0011822660098522167 
     0.00029556650246305416



### Location of the missing values

The positions of the missing data are evaluated by


```julia
mp = missingpos(mouse)
```




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
      [81   ,     5]  =  true
      [676  ,     5]  =  true
      [990  ,     5]  =  true
      ⋮
      [1791 , 10150]  =  true
      [1795 , 10150]  =  true
      [1846 , 10150]  =  true
      [1848 , 10150]  =  true
      [1851 , 10150]  =  true
      [1853 , 10150]  =  true
      [1860 , 10150]  =  true
      [1873 , 10150]  =  true
      [1886 , 10150]  =  true
      [1894 , 10150]  =  true
      [1897 , 10150]  =  true
      [1939 , 10150]  =  true




```julia
@btime(missingpos($mouse));
```

      35.144 ms (19273 allocations: 1.80 MiB)


So, for example, the number of missing data values in each column can be evaluated as


```julia
sum(mp, dims=1)
```




    1×10150 Array{Int64,2}:
     2  0  4  3  4  1  4  1  3  3  0  4  0  …  174  173  173  162  173  174  175



although it is faster, but somewhat more obscure, to use


```julia
view(counts(mouse, dims=1), 2:2, :)
```




    1×10150 view(::Array{Int64,2}, 2:2, :) with eltype Int64:
     2  0  4  3  4  1  4  1  3  3  0  4  0  …  174  173  173  162  173  174  175



## Genetic relationship matrix

`grm` function computes the empirical kinship matrix using either the classical genetic relationship matrix, `grm(A, model=:GRM)`, or the method of moment method, `grm(A, model=:MoM)`, or the robust method, `grm(A, model=:Robust)`. 

Classical genetic relation matrix


```julia
# grm(mouse, method=:MoM)
# grm(mouse, method=:Robust)
grm(mouse, method=:GRM)
```




    1940×1940 Array{Float64,2}:
      0.478301    -0.0331304    0.0135612    …  -0.0347737   -0.0129443 
     -0.0331304    0.422771    -0.0389227        0.0457987    0.00556832
      0.0135612   -0.0389227    0.509248        -0.0356689   -0.0608705 
      0.0198205    0.00728645  -0.00935362      -0.0302404   -0.0102152 
      0.056747    -0.0163418   -0.00495283      -0.0413347   -0.0415659 
     -0.0165628   -0.0191127   -0.0112181    …   0.0177118   -0.0193087 
      0.123771    -0.0404167    0.00442739       0.00880649  -0.0437565 
     -0.0628362    0.172552    -0.0728312        0.0640027   -0.0281429 
      0.0605018   -0.0260505    0.00398852      -0.00277754  -0.0607773 
      0.108886    -0.0204594   -0.00767711      -0.0210501    0.00343526
     -0.0142307    0.00270989  -0.0235504    …  -0.0223563   -0.028408  
     -0.0306022    0.197743    -0.00244269       0.0213998   -0.0478472 
     -0.0131463   -0.0226707    0.0223522       -0.037288     0.0493662 
      ⋮                                      ⋱                          
      0.0176725   -0.0165609    0.0378308        0.0238751   -0.0420143 
      0.0024949   -0.0411137    0.0154847       -0.0380656   -0.0650806 
      0.0952286    0.00894298  -0.0163446    …  -0.0202633   -0.0219594 
     -0.0309488   -0.0228342   -0.0478253       -0.014896     0.261623  
     -0.004804    -0.0375168   -0.0211418       -0.0172572    0.0359166 
      0.0076296    0.0481887   -0.0328968        0.0920425   -0.0292548 
      0.070045    -0.0302138    0.000647283      0.00892069  -0.00632566
      0.0378132   -6.59565e-5   0.00888932   …   0.00230815  -0.0291622 
     -0.00132837   0.00223654   0.0495928       -0.00936248   0.0299075 
      0.0640864   -0.0241218    0.00602283       0.00403413   0.00689551
     -0.0347737    0.0457987   -0.0356689        0.509228    -0.035215  
     -0.0129443    0.00556832  -0.0608705       -0.035215     0.552712  




```julia
@btime(grm($mouse, method=:GRM));
```

      466.965 ms (30 allocations: 28.95 MiB)


Using Float32 (single precision) potentially saves memory usage and computation time.


```julia
grm(mouse, method=:GRM, t=Float32)
```




    1940×1940 Array{Float32,2}:
      0.478301    -0.0331304    0.0135612    …  -0.0347737   -0.0129443 
     -0.0331304    0.422771    -0.0389227        0.0457987    0.00556833
      0.0135612   -0.0389227    0.509248        -0.0356689   -0.0608705 
      0.0198205    0.00728645  -0.00935361      -0.0302404   -0.0102152 
      0.056747    -0.0163418   -0.00495284      -0.0413347   -0.0415659 
     -0.0165628   -0.0191127   -0.0112181    …   0.0177117   -0.0193087 
      0.123771    -0.0404167    0.0044274        0.0088065   -0.0437565 
     -0.0628363    0.172552    -0.0728312        0.0640027   -0.0281429 
      0.0605018   -0.0260505    0.00398853      -0.00277754  -0.0607773 
      0.108886    -0.0204594   -0.00767711      -0.0210501    0.00343524
     -0.0142307    0.00270989  -0.0235504    …  -0.0223563   -0.028408  
     -0.0306022    0.197743    -0.00244268       0.0213998   -0.0478472 
     -0.0131464   -0.0226707    0.0223522       -0.037288     0.0493662 
      ⋮                                      ⋱                          
      0.0176725   -0.016561     0.0378308        0.0238751   -0.0420143 
      0.00249491  -0.0411137    0.0154847       -0.0380656   -0.0650806 
      0.0952286    0.00894298  -0.0163446    …  -0.0202633   -0.0219594 
     -0.0309488   -0.0228342   -0.0478253       -0.014896     0.261623  
     -0.00480401  -0.0375167   -0.0211418       -0.0172572    0.0359166 
      0.00762961   0.0481887   -0.0328968        0.0920425   -0.0292547 
      0.070045    -0.0302138    0.000647269      0.00892068  -0.00632566
      0.0378132   -6.59475e-5   0.00888932   …   0.00230815  -0.0291622 
     -0.00132838   0.00223653   0.0495928       -0.00936246   0.0299075 
      0.0640864   -0.0241219    0.00602283       0.00403413   0.00689551
     -0.0347737    0.0457987   -0.0356689        0.509228    -0.035215  
     -0.0129443    0.00556833  -0.0608705       -0.035215     0.552712  




```julia
@btime(grm($mouse, method=:GRM, t=Float32));
```

      276.116 ms (31 allocations: 14.60 MiB)


By default, `grm` exlcude SNPs with minor allele frequency below 0.01. This can be changed by the keyword argument `minmaf`.


```julia
# compute GRM excluding SNPs with MAF≤0.05 
grm(mouse, minmaf=0.05)
```




    1940×1940 Array{Float64,2}:
      0.478556    -0.0331783    0.013541     …  -0.0348225   -0.0129761 
     -0.0331783    0.422993    -0.0389741        0.0457975    0.00554753
      0.013541    -0.0389741    0.50952         -0.0357183   -0.0609305 
      0.0203209    0.00777944  -0.00887047      -0.0297696   -0.00972836
      0.0567523   -0.0163798   -0.00498406      -0.0413874   -0.0416146 
     -0.0166009   -0.0191523   -0.0112531    …   0.0176939   -0.0193442 
      0.123816    -0.0404689    0.00440171       0.0087834   -0.0438065 
     -0.0629017    0.172626    -0.0729026        0.0640123   -0.0281836 
      0.0605093   -0.0260942    0.00396257      -0.00280748  -0.0608373 
      0.108922    -0.0204998   -0.00770996      -0.0210909    0.00341321
     -0.0142674    0.00268319  -0.0235927    …  -0.0223978   -0.0284489 
     -0.0306486    0.197832    -0.00247243       0.0213842   -0.0478996 
     -0.0131824   -0.0227124    0.0223371       -0.0373384    0.0493713 
      ⋮                                      ⋱                          
      0.0176546   -0.016599     0.0378249        0.0238609   -0.0420633 
      0.00246808  -0.0411663    0.0154656       -0.0381165   -0.0651432 
      0.0952566    0.00891997  -0.0163826    …  -0.0203036   -0.0219965 
     -0.0309912   -0.0228718   -0.0478777       -0.0149289    0.261754  
     -0.00483514  -0.0375673   -0.0211827       -0.0172957    0.0359138 
      0.00770862   0.0482917   -0.0328417        0.0921714   -0.0292961 
      0.0700582   -0.03026      0.000619365      0.00889767  -0.00635348
      0.0378313   -7.02155e-5   0.00889036   …   0.0023053   -0.0291795 
     -0.00133338   0.00223364   0.0496179       -0.00937223   0.0299252 
      0.0641201   -0.0241403    0.00602217       0.0040323    0.00689958
     -0.0348225    0.0457975   -0.0357183        0.509501    -0.0352599 
     -0.0129761    0.00554753  -0.0609305       -0.0352599    0.553015  



To specify specific SNPs for calculating empirical kinship, use the `cinds` keyword (default is `nothing`). When `cinds` is specified, `minmaf` is ignored.


```julia
# GRM using every other SNP
grm(mouse, cinds=1:2:size(mouse, 2))
```




    1940×1940 Array{Float64,2}:
      0.477       -0.0307774     0.0118026   …  -0.0320301    -0.0125113 
     -0.0307774    0.425085     -0.0367459       0.0480442     0.00519065
      0.0118026   -0.0367459     0.505038       -0.0385129    -0.0631557 
      0.0166017    0.00614789   -0.00919695     -0.0399744    -0.0104884 
      0.05724     -0.0122148    -0.00543377     -0.0395663    -0.0372998 
     -0.0193129   -0.0224378    -0.009277    …   0.0153785    -0.0220184 
      0.12194     -0.0410682     0.00274307      0.00796748   -0.0441578 
     -0.0624031    0.173985     -0.0724784       0.0663191    -0.0294243 
      0.0627626   -0.0288615     0.00265615     -0.00449877   -0.0579702 
      0.110878    -0.0232715    -0.00881604     -0.021272      0.00169016
     -0.00800735  -0.00149824   -0.019791    …  -0.024124     -0.0289397 
     -0.0272944    0.19894      -0.00534771      0.0209384    -0.0511051 
     -0.011388    -0.0281003     0.0273853      -0.0360047     0.0459359 
      ⋮                                      ⋱                           
      0.0169431   -0.0136989     0.0340794       0.0272811    -0.041189  
      0.00201325  -0.0426611     0.0124353      -0.0387982    -0.0656181 
      0.097587     0.0058123    -0.0160698   …  -0.021457     -0.023226  
     -0.0342014   -0.0211246    -0.0490112      -0.0129575     0.256552  
     -0.00324255  -0.0423482    -0.0192699      -0.0149015     0.0339388 
      0.00575353   0.0464237    -0.0294694       0.0924759    -0.0275451 
      0.0748725   -0.0258461    -0.00141068      0.0115232    -0.00486589
      0.0386555    0.000612169   0.00959997  …  -0.000357284  -0.0334687 
     -0.00343056   0.0120673     0.0455375      -0.0103798     0.0336959 
      0.0656909   -0.0193469     0.00600815      0.00188545    0.00726181
     -0.0320301    0.0480442    -0.0385129       0.513285     -0.0317963 
     -0.0125113    0.00519065   -0.0631557      -0.0317963     0.54471   



## Fitering

Before GWAS, we often need to filter SNPs and/or samples according to genotyping success rates, minor allele frequencies, and Hardy-Weinberg Equilibrium test. This can be achieved by the `filter` function.

```@docs
SnpArrays.filter
```

By default, it outputs row and column index vectors such that sample-wise and SNP-wise genotyping success rate are at least 0.98 and minor allele frequencies are at least 0.01. User can opt to filter according to Hardy-Weinberg test by setting the minumum p-value `min_hwe_pval`.


```julia
rowmask, colmask =  SnpArrays.filter(mouse)
```




    (Bool[true, true, true, true, true, true, true, true, true, true  …  true, true, true, true, true, true, true, true, true, true], Bool[true, true, true, true, true, true, true, true, true, true  …  false, false, false, false, false, false, false, false, false, false])




```julia
count(rowmask), count(colmask)
```




    (1930, 10072)




```julia
@btime(SnpArrays.filter($mouse, min_success_rate_per_row=0.999, min_success_rate_per_col=0.999));
```

      85.756 ms (7 allocations: 17.17 KiB)


One may use the `rowmask` and `colmask` to filter and save filtering result as Plink files.
```julia
SnpArrays.filter(SnpArrays.datadir("mouse"), rowmask, colmask)
```

## Filter Plink files

Filter a set of Plink files according to row indices and column indices. By result, filtered Plink files are saved as `srcname.filtered.bed`, `srcname.filtered.fam`, and `srcname.filtered.bim`, where `srcname` is the source Plink file name. You can also specify destimation file name using keyword `des`.


```julia
SnpArrays.filter(SnpArrays.datadir("mouse"), 1:5, 1:5)
```




    5×5 SnpArray:
     0x02  0x02  0x02  0x02  0x03
     0x02  0x02  0x03  0x02  0x02
     0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02
     0x03  0x03  0x03  0x03  0x03




```julia
# clean up
rm(SnpArrays.datadir("mouse.filtered.bed"), force=true)
rm(SnpArrays.datadir("mouse.filtered.fam"), force=true)
rm(SnpArrays.datadir("mouse.filtered.bim"), force=true)
```

Filter a set of Plink files according to logical vectors.


```julia
SnpArrays.filter(SnpArrays.datadir("mouse"), rowmask, colmask)
```




    1930×10072 SnpArray:
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x02  0x03  0x02  0x02  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x00  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x00  0x00  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02     0x03  0x03  0x03  0x00  0x00  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x02  0x03  0x02  0x02  0x02  0x03
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x00  0x00  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x00  0x03  0x00  0x00
     0x02  0x02  0x03  0x02  0x02  0x02     0x00  0x03  0x03  0x00  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x00  0x03  0x00  0x00  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x02  0x03  0x02  0x02  0x02  0x03
     0x03  0x03  0x03  0x03  0x03  0x03  …  0x03  0x03  0x00  0x03  0x00  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x00  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x00  0x00  0x03  0x03  0x03
        ⋮                             ⋮  ⋱                             ⋮      
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x00  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x00  0x00  0x00
     0x02  0x02  0x02  0x02  0x02  0x02  …  0x03  0x03  0x03  0x00  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x00  0x03  0x00  0x00  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x00  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x00  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x02  0x03  0x03  0x00  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x03  0x00  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x00  0x00  0x03  0x00  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x00  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02     0x03  0x00  0x00  0x03  0x03  0x03
     0x00  0x00  0x00  0x00  0x03  0x00     0x03  0x03  0x00  0x03  0x00  0x03




```julia
readdir(glob"mouse.filtered.*", datapath)
```




    3-element Array{String,1}:
     "/Users/huazhou/.julia/dev/SnpArrays/data/mouse.filtered.bed"
     "/Users/huazhou/.julia/dev/SnpArrays/data/mouse.filtered.bim"
     "/Users/huazhou/.julia/dev/SnpArrays/data/mouse.filtered.fam"




```julia
# clean up
rm(SnpArrays.datadir("mouse.filtered.bed"), force=true)
rm(SnpArrays.datadir("mouse.filtered.fam"), force=true)
rm(SnpArrays.datadir("mouse.filtered.bim"), force=true)
```

# SnpBitMatrix

In some applications we want to perform linear algebra using SnpArray directly without expanding it to numeric matrix. This is achieved by the `SnpBitMatrix` type. The implementation assumes:

1. the SnpArray does not have missing genotypes, and
2. the matrix corresponding to SnpArray is the matrix of A2 allele counts.

## Constructor

First let's load a data set without missing genotypes.


```julia
const EUR = SnpArray(SnpArrays.datadir("EUR_subset.bed"))
```




    379×54051 SnpArray:
     0x03  0x03  0x03  0x02  0x02  0x03  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x02  0x03  0x02  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x02  0x02  0x02  0x03  0x03  0x02
     0x03  0x03  0x03  0x00  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x00  0x03  0x03     0x02  0x02  0x02  0x03  0x03  0x03
     0x02  0x03  0x03  0x03  0x03  0x03  …  0x03  0x03  0x03  0x03  0x03  0x02
     0x02  0x03  0x03  0x02  0x02  0x03     0x03  0x03  0x02  0x02  0x03  0x03
     0x02  0x03  0x03  0x03  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x00  0x02  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x03  0x03  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x02  0x03  0x03  …  0x03  0x03  0x02  0x02  0x03  0x03
     0x03  0x03  0x03  0x02  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x02
     0x03  0x02  0x03  0x02  0x02  0x03     0x03  0x03  0x03  0x03  0x03  0x03
        ⋮                             ⋮  ⋱     ⋮                             ⋮
     0x03  0x03  0x03  0x00  0x02  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x02  0x02  0x03     0x02  0x02  0x02  0x03  0x02  0x03
     0x03  0x03  0x03  0x02  0x02  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x03  0x03  0x02  0x03  0x03  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x00  0x00  0x03     0x02  0x02  0x02  0x03  0x03  0x03
     0x02  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x02  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x02  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x03  0x03  0x03  0x03  0x03  …  0x03  0x03  0x02  0x02  0x03  0x03
     0x03  0x03  0x03  0x00  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x03  0x03  0x02  0x00  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x02  0x02  0x03     0x03  0x03  0x03  0x03  0x03  0x03



To instantiate a SnpBitMatrix based on SnpArray,


```julia
const EURbm = SnpBitMatrix{Float64}(EUR, model=ADDITIVE_MODEL, center=true, scale=true);
```

The constructor shares the same keyword arguments as the `convert` or `copyto!` functions. The type parameter, `Float64` in this example, indicates the SnpBitMatrix acts like a Float64 matrix.

The memory usage of the SnpBitMatrix should be similar to the SnpArray, or equivalently bed file size, if `model=ADDITIVE_MODEL`, or half of that of SnpArray if `model=DOMINANT_MODEL` or `model=RECESSIVE_MODEL`.


```julia
Base.summarysize(EUR), Base.summarysize(EURbm)
```




    (6876757, 6421960)



## Linear algebra

A SnpBitMatrix acts similar to a regular matrix and responds to `size`, `eltype`, and SnpBitMatrix-vector multiplications.


```julia
@show size(EURbm)
@show eltype(EURbm)
@show typeof(EURbm) <: AbstractMatrix;
```

    size(EURbm) = (379, 54051)
    eltype(EURbm) = Float64
    typeof(EURbm) <: AbstractMatrix = true


SnpBitMatrix-vector multiplication is mathematically equivalent to the corresponding Float matrix contained from `convert` or `copyto!` a SnpArray.


```julia
using LinearAlgebra
v1 = randn(size(EUR, 1))
v2 = randn(size(EUR, 2))
A = convert(Matrix{Float64}, EUR, model=ADDITIVE_MODEL, center=true, scale=true)
norm(EURbm * v2 -  A * v2)
```




    6.484489785912384e-11




```julia
norm(EURbm' * v1 - A' * v1)
```




    8.354707200357845e-12



In this example, the Float64 matrix fits into memory so the SnpBitMatrix-vector multiplication is much slower than Matrix{Float64}-vector multiplication (highly optimized BLAS).


```julia
out1 = Vector{Float64}(undef, size(EUR, 1))
out2 = Vector{Float64}(undef, size(EUR, 2))
@btime(mul!($out1, $EURbm, $v2));
```

      83.763 ms (0 allocations: 0 bytes)



```julia
@btime(mul!($out1, $A, $v2));
```

      7.727 ms (0 allocations: 0 bytes)



```julia
@btime(mul!($out2, $transpose($EURbm), $v1));
```

      78.741 ms (1 allocation: 16 bytes)



```julia
@btime(mul!($out2, $transpose($A), $v1));
```

      6.715 ms (0 allocations: 0 bytes)


In another test example with ~1GB bed file, SnpBitMatrix-vector multiplication is about 3-5 folder faster than the corresponding Matrix{Float64}-vector multiplication, because the Matrix{Float64} matrix cannot fit into the memory.

`SnpBitMatrix` can be created from a subarray of SnpArray.


```julia
EURsub = @view EUR[1:2:100, 1:2:100]
EURsubbm = SnpBitMatrix{Float64}(EURsub, model=ADDITIVE_MODEL, center=true, scale=true);
```


```julia
Base.summarysize(EURsubbm)
```




    2600




```julia
@show size(EURsubbm)
@show eltype(EURsubbm)
@show typeof(EURsubbm) <: AbstractMatrix;
```

    size(EURsubbm) = (50, 50)
    eltype(EURsubbm) = Float64
    typeof(EURsubbm) <: AbstractMatrix = true



```julia
using LinearAlgebra
v1 = randn(size(EURsub, 1))
v2 = randn(size(EURsub, 2))
A = convert(Matrix{Float64}, EURsub, model=ADDITIVE_MODEL, center=true, scale=true)
norm(EURsubbm * v2 -  A * v2)
```




    3.692455572052395e-14




```julia
norm(EURsubbm' * v1 - A' * v1)
```




    1.460805819955132e-13



# SnpData

We can create a `SnpData`, which has a `SnpArray` with information on SNP and subject appended.

## Constructor


```julia
EUR_data = SnpData(SnpArrays.datadir("EUR_subset"))
```




    SnpData(people: 379, snps: 54051,
    snp_info: 
    │ Row │ chromosome │ snpid       │ genetic_distance │ position │ allele1      │ allele2      │
    │     │ String     │ String      │ Float64          │ Int64    │ Categorical… │ Categorical… │
    ├─────┼────────────┼─────────────┼──────────────────┼──────────┼──────────────┼──────────────┤
    │ 1   │ 17         │ rs34151105  │ 0.0              │ 1665     │ T            │ C            │
    │ 2   │ 17         │ rs143500173 │ 0.0              │ 2748     │ T            │ A            │
    │ 3   │ 17         │ rs113560219 │ 0.0              │ 4702     │ T            │ C            │
    │ 4   │ 17         │ rs1882989   │ 5.6e-5           │ 15222    │ G            │ A            │
    │ 5   │ 17         │ rs8069133   │ 0.000499         │ 32311    │ G            │ A            │
    │ 6   │ 17         │ rs112221137 │ 0.000605         │ 36405    │ G            │ T            │
    …,
    person_info: 
    │ Row │ fid       │ iid       │ father    │ mother    │ sex       │ phenotype │
    │     │ Abstract… │ Abstract… │ Abstract… │ Abstract… │ Abstract… │ Abstract… │
    ├─────┼───────────┼───────────┼───────────┼───────────┼───────────┼───────────┤
    │ 1   │ 1         │ HG00096   │ 0         │ 0         │ 1         │ 1         │
    │ 2   │ 2         │ HG00097   │ 0         │ 0         │ 2         │ 1         │
    │ 3   │ 3         │ HG00099   │ 0         │ 0         │ 2         │ 1         │
    │ 4   │ 4         │ HG00100   │ 0         │ 0         │ 2         │ 1         │
    │ 5   │ 5         │ HG00101   │ 0         │ 0         │ 1         │ 1         │
    │ 6   │ 6         │ HG00102   │ 0         │ 0         │ 2         │ 1         │
    …,
    src: /Users/huazhou/.julia/dev/SnpArrays/src/../data/EUR_subset
    )



## Filter

We can filter SnpData by functions `f_person` and `f_snp`. `f_person` applies to the field `person_info` and selects persons (rows) for which `f_person` is `true`.`f_snp` applies to the field `snp_info` and selects snps (columns) for which `f_snp` is `true`. The first argument can be either a `SnpData` or an `AbstractString`.


```julia
SnpArrays.filter(EUR_data; des="tmp.filter.chr.17", f_snp = x -> x[:chromosome]=="17")
```




    SnpData(people: 379, snps: 11041,
    snp_info: 
    │ Row │ chromosome │ snpid       │ genetic_distance │ position │ allele1      │ allele2      │
    │     │ String     │ String      │ Float64          │ Int64    │ Categorical… │ Categorical… │
    ├─────┼────────────┼─────────────┼──────────────────┼──────────┼──────────────┼──────────────┤
    │ 1   │ 17         │ rs34151105  │ 0.0              │ 1665     │ T            │ C            │
    │ 2   │ 17         │ rs143500173 │ 0.0              │ 2748     │ T            │ A            │
    │ 3   │ 17         │ rs113560219 │ 0.0              │ 4702     │ T            │ C            │
    │ 4   │ 17         │ rs1882989   │ 5.6e-5           │ 15222    │ G            │ A            │
    │ 5   │ 17         │ rs8069133   │ 0.000499         │ 32311    │ G            │ A            │
    │ 6   │ 17         │ rs112221137 │ 0.000605         │ 36405    │ G            │ T            │
    …,
    person_info: 
    │ Row │ fid       │ iid       │ father    │ mother    │ sex       │ phenotype │
    │     │ Abstract… │ Abstract… │ Abstract… │ Abstract… │ Abstract… │ Abstract… │
    ├─────┼───────────┼───────────┼───────────┼───────────┼───────────┼───────────┤
    │ 1   │ 1         │ HG00096   │ 0         │ 0         │ 1         │ 1         │
    │ 2   │ 2         │ HG00097   │ 0         │ 0         │ 2         │ 1         │
    │ 3   │ 3         │ HG00099   │ 0         │ 0         │ 2         │ 1         │
    │ 4   │ 4         │ HG00100   │ 0         │ 0         │ 2         │ 1         │
    │ 5   │ 5         │ HG00101   │ 0         │ 0         │ 1         │ 1         │
    │ 6   │ 6         │ HG00102   │ 0         │ 0         │ 2         │ 1         │
    …,
    src: tmp.filter.chr.17
    )




```julia
SnpArrays.filter(SnpArrays.datadir("EUR_subset"); des="tmp.filter.chr.17", f_snp = x -> x[:chromosome]=="17")
```




    SnpData(people: 379, snps: 11041,
    snp_info: 
    │ Row │ chromosome │ snpid       │ genetic_distance │ position │ allele1      │ allele2      │
    │     │ String     │ String      │ Float64          │ Int64    │ Categorical… │ Categorical… │
    ├─────┼────────────┼─────────────┼──────────────────┼──────────┼──────────────┼──────────────┤
    │ 1   │ 17         │ rs34151105  │ 0.0              │ 1665     │ T            │ C            │
    │ 2   │ 17         │ rs143500173 │ 0.0              │ 2748     │ T            │ A            │
    │ 3   │ 17         │ rs113560219 │ 0.0              │ 4702     │ T            │ C            │
    │ 4   │ 17         │ rs1882989   │ 5.6e-5           │ 15222    │ G            │ A            │
    │ 5   │ 17         │ rs8069133   │ 0.000499         │ 32311    │ G            │ A            │
    │ 6   │ 17         │ rs112221137 │ 0.000605         │ 36405    │ G            │ T            │
    …,
    person_info: 
    │ Row │ fid       │ iid       │ father    │ mother    │ sex       │ phenotype │
    │     │ Abstract… │ Abstract… │ Abstract… │ Abstract… │ Abstract… │ Abstract… │
    ├─────┼───────────┼───────────┼───────────┼───────────┼───────────┼───────────┤
    │ 1   │ 1         │ HG00096   │ 0         │ 0         │ 1         │ 1         │
    │ 2   │ 2         │ HG00097   │ 0         │ 0         │ 2         │ 1         │
    │ 3   │ 3         │ HG00099   │ 0         │ 0         │ 2         │ 1         │
    │ 4   │ 4         │ HG00100   │ 0         │ 0         │ 2         │ 1         │
    │ 5   │ 5         │ HG00101   │ 0         │ 0         │ 1         │ 1         │
    │ 6   │ 6         │ HG00102   │ 0         │ 0         │ 2         │ 1         │
    …,
    src: tmp.filter.chr.17
    )




```julia
SnpArrays.filter(EUR_data; des="tmp.filter.sex.male", f_person = x -> x[:sex] == "1")
```




    SnpData(people: 178, snps: 54051,
    snp_info: 
    │ Row │ chromosome │ snpid       │ genetic_distance │ position │ allele1      │ allele2      │
    │     │ String     │ String      │ Float64          │ Int64    │ Categorical… │ Categorical… │
    ├─────┼────────────┼─────────────┼──────────────────┼──────────┼──────────────┼──────────────┤
    │ 1   │ 17         │ rs34151105  │ 0.0              │ 1665     │ T            │ C            │
    │ 2   │ 17         │ rs143500173 │ 0.0              │ 2748     │ T            │ A            │
    │ 3   │ 17         │ rs113560219 │ 0.0              │ 4702     │ T            │ C            │
    │ 4   │ 17         │ rs1882989   │ 5.6e-5           │ 15222    │ G            │ A            │
    │ 5   │ 17         │ rs8069133   │ 0.000499         │ 32311    │ G            │ A            │
    │ 6   │ 17         │ rs112221137 │ 0.000605         │ 36405    │ G            │ T            │
    …,
    person_info: 
    │ Row │ fid       │ iid       │ father    │ mother    │ sex       │ phenotype │
    │     │ Abstract… │ Abstract… │ Abstract… │ Abstract… │ Abstract… │ Abstract… │
    ├─────┼───────────┼───────────┼───────────┼───────────┼───────────┼───────────┤
    │ 1   │ 1         │ HG00096   │ 0         │ 0         │ 1         │ 1         │
    │ 2   │ 5         │ HG00101   │ 0         │ 0         │ 1         │ 1         │
    │ 3   │ 7         │ HG00103   │ 0         │ 0         │ 1         │ 1         │
    │ 4   │ 10        │ HG00108   │ 0         │ 0         │ 1         │ 1         │
    │ 5   │ 11        │ HG00109   │ 0         │ 0         │ 1         │ 1         │
    │ 6   │ 14        │ HG00112   │ 0         │ 0         │ 1         │ 1         │
    …,
    src: tmp.filter.sex.male
    )



Both `f_person` and `f_snp` can be used at the same time.


```julia
SnpArrays.filter(EUR_data; des="tmp.filter.chr.17.sex.male", f_person = x -> x[:sex] == "1", f_snp = x -> x[:chromosome] == "17")
```




    SnpData(people: 178, snps: 11041,
    snp_info: 
    │ Row │ chromosome │ snpid       │ genetic_distance │ position │ allele1      │ allele2      │
    │     │ String     │ String      │ Float64          │ Int64    │ Categorical… │ Categorical… │
    ├─────┼────────────┼─────────────┼──────────────────┼──────────┼──────────────┼──────────────┤
    │ 1   │ 17         │ rs34151105  │ 0.0              │ 1665     │ T            │ C            │
    │ 2   │ 17         │ rs143500173 │ 0.0              │ 2748     │ T            │ A            │
    │ 3   │ 17         │ rs113560219 │ 0.0              │ 4702     │ T            │ C            │
    │ 4   │ 17         │ rs1882989   │ 5.6e-5           │ 15222    │ G            │ A            │
    │ 5   │ 17         │ rs8069133   │ 0.000499         │ 32311    │ G            │ A            │
    │ 6   │ 17         │ rs112221137 │ 0.000605         │ 36405    │ G            │ T            │
    …,
    person_info: 
    │ Row │ fid       │ iid       │ father    │ mother    │ sex       │ phenotype │
    │     │ Abstract… │ Abstract… │ Abstract… │ Abstract… │ Abstract… │ Abstract… │
    ├─────┼───────────┼───────────┼───────────┼───────────┼───────────┼───────────┤
    │ 1   │ 1         │ HG00096   │ 0         │ 0         │ 1         │ 1         │
    │ 2   │ 5         │ HG00101   │ 0         │ 0         │ 1         │ 1         │
    │ 3   │ 7         │ HG00103   │ 0         │ 0         │ 1         │ 1         │
    │ 4   │ 10        │ HG00108   │ 0         │ 0         │ 1         │ 1         │
    │ 5   │ 11        │ HG00109   │ 0         │ 0         │ 1         │ 1         │
    │ 6   │ 14        │ HG00112   │ 0         │ 0         │ 1         │ 1         │
    …,
    src: tmp.filter.chr.17.sex.male
    )



## Split

We can split `SnpData` by SNP's choromosomes or each person's sex or phenotype using `split_plink`. Again, the first argument can be an `SnpData` or an `AbstractString`.


```julia
splitted = SnpArrays.split_plink(SnpArrays.datadir("EUR_subset"), :chromosome; prefix="tmp.split.chr.")
```




    Dict{AbstractString,SnpData} with 6 entries:
      "21" => SnpData(people: 379, snps: 5813,…
      "17" => SnpData(people: 379, snps: 11041,…
      "19" => SnpData(people: 379, snps: 9690,…
      "20" => SnpData(people: 379, snps: 9327,…
      "22" => SnpData(people: 379, snps: 5938,…
      "18" => SnpData(people: 379, snps: 12242,…



Let's take a SnpArray for chromosome 17.


```julia
piece = splitted["17"]
```




    SnpData(people: 379, snps: 11041,
    snp_info: 
    │ Row │ chromosome │ snpid       │ genetic_distance │ position │ allele1      │ allele2      │
    │     │ String     │ String      │ Float64          │ Int64    │ Categorical… │ Categorical… │
    ├─────┼────────────┼─────────────┼──────────────────┼──────────┼──────────────┼──────────────┤
    │ 1   │ 17         │ rs34151105  │ 0.0              │ 1665     │ T            │ C            │
    │ 2   │ 17         │ rs143500173 │ 0.0              │ 2748     │ T            │ A            │
    │ 3   │ 17         │ rs113560219 │ 0.0              │ 4702     │ T            │ C            │
    │ 4   │ 17         │ rs1882989   │ 5.6e-5           │ 15222    │ G            │ A            │
    │ 5   │ 17         │ rs8069133   │ 0.000499         │ 32311    │ G            │ A            │
    │ 6   │ 17         │ rs112221137 │ 0.000605         │ 36405    │ G            │ T            │
    …,
    person_info: 
    │ Row │ fid       │ iid       │ father    │ mother    │ sex       │ phenotype │
    │     │ Abstract… │ Abstract… │ Abstract… │ Abstract… │ Abstract… │ Abstract… │
    ├─────┼───────────┼───────────┼───────────┼───────────┼───────────┼───────────┤
    │ 1   │ 1         │ HG00096   │ 0         │ 0         │ 1         │ 1         │
    │ 2   │ 2         │ HG00097   │ 0         │ 0         │ 2         │ 1         │
    │ 3   │ 3         │ HG00099   │ 0         │ 0         │ 2         │ 1         │
    │ 4   │ 4         │ HG00100   │ 0         │ 0         │ 2         │ 1         │
    │ 5   │ 5         │ HG00101   │ 0         │ 0         │ 1         │ 1         │
    │ 6   │ 6         │ HG00102   │ 0         │ 0         │ 2         │ 1         │
    …,
    src: tmp.split.chr.17
    )




```julia
@assert all(piece.snp_info[:chromosome].== "17")
```


```julia
splitted_sex = SnpArrays.split_plink(EUR_data, :sex; prefix="tmp.split.sex.")
```




    Dict{AbstractString,SnpData} with 2 entries:
      "1" => SnpData(people: 178, snps: 54051,…
      "2" => SnpData(people: 201, snps: 54051,…



## Merge

We can merge the splitted dictionary back into one SnpData using `merge_plink`.


```julia
merged = SnpArrays.merge_plink("tmp.merged", splitted) # write_plink is included here
```




    SnpData(people: 379, snps: 54051,
    snp_info: 
    │ Row │ chromosome │ snpid       │ genetic_distance │ position │ allele1      │ allele2      │
    │     │ String     │ String      │ Float64          │ Int64    │ Categorical… │ Categorical… │
    ├─────┼────────────┼─────────────┼──────────────────┼──────────┼──────────────┼──────────────┤
    │ 1   │ 17         │ rs34151105  │ 0.0              │ 1665     │ T            │ C            │
    │ 2   │ 17         │ rs143500173 │ 0.0              │ 2748     │ T            │ A            │
    │ 3   │ 17         │ rs113560219 │ 0.0              │ 4702     │ T            │ C            │
    │ 4   │ 17         │ rs1882989   │ 5.6e-5           │ 15222    │ G            │ A            │
    │ 5   │ 17         │ rs8069133   │ 0.000499         │ 32311    │ G            │ A            │
    │ 6   │ 17         │ rs112221137 │ 0.000605         │ 36405    │ G            │ T            │
    …,
    person_info: 
    │ Row │ fid       │ iid       │ father    │ mother    │ sex       │ phenotype │
    │     │ Abstract… │ Abstract… │ Abstract… │ Abstract… │ Abstract… │ Abstract… │
    ├─────┼───────────┼───────────┼───────────┼───────────┼───────────┼───────────┤
    │ 1   │ 1         │ HG00096   │ 0         │ 0         │ 1         │ 1         │
    │ 2   │ 2         │ HG00097   │ 0         │ 0         │ 2         │ 1         │
    │ 3   │ 3         │ HG00099   │ 0         │ 0         │ 2         │ 1         │
    │ 4   │ 4         │ HG00100   │ 0         │ 0         │ 2         │ 1         │
    │ 5   │ 5         │ HG00101   │ 0         │ 0         │ 1         │ 1         │
    │ 6   │ 6         │ HG00102   │ 0         │ 0         │ 2         │ 1         │
    …,
    src: tmp.merged
    )



You can also merge the plink formatted files based on their common prefix.


```julia
merged_from_splitted_files = merge_plink("tmp.split.chr"; des = "tmp.merged.2")
```




    SnpData(people: 379, snps: 54051,
    snp_info: 
    │ Row │ chromosome │ snpid       │ genetic_distance │ position │ allele1      │ allele2      │
    │     │ String     │ String      │ Float64          │ Int64    │ Categorical… │ Categorical… │
    ├─────┼────────────┼─────────────┼──────────────────┼──────────┼──────────────┼──────────────┤
    │ 1   │ 17         │ rs34151105  │ 0.0              │ 1665     │ T            │ C            │
    │ 2   │ 17         │ rs143500173 │ 0.0              │ 2748     │ T            │ A            │
    │ 3   │ 17         │ rs113560219 │ 0.0              │ 4702     │ T            │ C            │
    │ 4   │ 17         │ rs1882989   │ 5.6e-5           │ 15222    │ G            │ A            │
    │ 5   │ 17         │ rs8069133   │ 0.000499         │ 32311    │ G            │ A            │
    │ 6   │ 17         │ rs112221137 │ 0.000605         │ 36405    │ G            │ T            │
    …,
    person_info: 
    │ Row │ fid       │ iid       │ father    │ mother    │ sex       │ phenotype │
    │     │ Abstract… │ Abstract… │ Abstract… │ Abstract… │ Abstract… │ Abstract… │
    ├─────┼───────────┼───────────┼───────────┼───────────┼───────────┼───────────┤
    │ 1   │ 1         │ HG00096   │ 0         │ 0         │ 1         │ 1         │
    │ 2   │ 2         │ HG00097   │ 0         │ 0         │ 2         │ 1         │
    │ 3   │ 3         │ HG00099   │ 0         │ 0         │ 2         │ 1         │
    │ 4   │ 4         │ HG00100   │ 0         │ 0         │ 2         │ 1         │
    │ 5   │ 5         │ HG00101   │ 0         │ 0         │ 1         │ 1         │
    │ 6   │ 6         │ HG00102   │ 0         │ 0         │ 2         │ 1         │
    …,
    src: tmp.merged.2
    )




```julia
# clean up
for ft in ["bim", "fam", "bed"]
    rm("tmp.filter.chr.17." * ft, force=true)
    rm("tmp.filter.sex.male." * ft, force=true)
    rm("tmp.filter.chr.17.sex.male." * ft, force=true)
    for k in keys(splitted)
        rm("tmp.split.chr.$(k)." * ft, force=true)
    end
    for k in keys(splitted_sex)
        rm("tmp.split.sex.$(k)." * ft, force=true)
    end
    rm("tmp.merged." * ft, force=true)
    rm("tmp.merged.2." * ft, force=true)
end
```
