
# SnpArray

`SnpArray` is an array of `Tuple{Bool,Bool}` and adopts the same coding as the [Plink binary format](http://zzz.bwh.harvard.edu/plink/binary.shtml). If `A1` and `A2` are the two alleles, the coding rule is  

| Genotype | SnpArray |  
|:---:|:---:|  
| A1,A1 | (false,false) |  
| A1,A2 | (false,true) |  
| A2,A2 | (true,true) |  
| missing | (true,false) |  

The code `(true,false)` is reserved for missing genotype. Otherwise, the bit `true` represents one copy of allele `A2`. In a two-dimensional `SnpArray`, each row is a person and each column is a SNP.

For complete genotype data, for example, after imputation, consider using the [HaplotypeArray](@ref) type.

## Constructor

There are various ways to initialize a `SnpArray`.  

* `SnpArray` can be initialized from [Plink binary files](http://zzz.bwh.harvard.edu/plink/binary.shtml), say the sample data set `hapmap3` in package `docs` folder:


```julia
;ls -l $(Pkg.dir("SnpArrays") * "/docs/hapmap3.*")
```

    -rw-r--r--  1 huazhou  staff  1128171 Jun 19  2017 /Users/huazhou/.julia/v0.6/SnpArrays/docs/hapmap3.bed
    -rw-r--r--  1 huazhou  staff   388672 Jun 19  2017 /Users/huazhou/.julia/v0.6/SnpArrays/docs/hapmap3.bim
    -rw-r--r--  1 huazhou  staff     7136 Jun 19  2017 /Users/huazhou/.julia/v0.6/SnpArrays/docs/hapmap3.fam
    -rw-r--r--  1 huazhou  staff   332960 Jun 19  2017 /Users/huazhou/.julia/v0.6/SnpArrays/docs/hapmap3.map



```julia
using SnpArrays
hapmap = SnpArray(Pkg.dir("SnpArrays") * "/docs/hapmap3")
```

    [1m[36mINFO: [39m[22m[36mv1.0 BED file detected
    [39m




    324×13928 SnpArrays.SnpArray{2}:
     (true, true)  (true, true)   (false, false)  …  (true, true)   (true, true)
     (true, true)  (false, true)  (false, true)      (false, true)  (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (false, true)  (true, true)    …  (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)     (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)   …  (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)     (true, true)   (true, true)
     ⋮                                            ⋱                             
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)     (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)     (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)  …  (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)  …  (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)     (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)     (true, true)   (true, true)



By default, the constructor figures out the number of individuals and SNPs from the `.bim` and `.fam` files.


```julia
# rows are people; columns are SNPs
people, snps = size(hapmap)
```




    (324, 13928)



Alternatively, users can supply keyword arguments `people` and `snps` to the constructor. In this case only the `.bed` file needs to be present.


```julia
hapmap = SnpArray(Pkg.dir("SnpArrays") * "/docs/hapmap3"; people = 324, snps = 13928)
```

    [1m[36mINFO: [39m[22m[36mv1.0 BED file detected
    [39m




    324×13928 SnpArrays.SnpArray{2}:
     (true, true)  (true, true)   (false, false)  …  (true, true)   (true, true)
     (true, true)  (false, true)  (false, true)      (false, true)  (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (false, true)  (true, true)    …  (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)     (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)   …  (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)     (true, true)   (true, true)
     ⋮                                            ⋱                             
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)     (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)     (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)  …  (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)  …  (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)     (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)     (true, true)   (true, true)



Internally `SnpArray` stores data as `BitArray`s and consumes approximately the same amount of memory as the Plink `bed` file size.


```julia
# memory usage, bed file size
Base.summarysize(hapmap), filesize("hapmap3.bed")
```




    (1128256, 1128171)



* `SnpArray` can be initialized from a matrix of A1 allele counts.


```julia
SnpArray(rand(0:2, 5, 3))
```




    5×3 SnpArrays.SnpArray{2}:
     (true, true)    (true, true)   (false, false)
     (false, false)  (true, true)   (false, true) 
     (false, false)  (false, true)  (false, true) 
     (false, false)  (false, true)  (true, true)  
     (true, true)    (true, true)   (true, true)  



* `SnpArray(m, n)` generates an m by n `SnpArray` of all A1 alleles.


```julia
s = SnpArray(5, 3)
```




    5×3 SnpArrays.SnpArray{2}:
     (false, false)  (false, false)  (false, false)
     (false, false)  (false, false)  (false, false)
     (false, false)  (false, false)  (false, false)
     (false, false)  (false, false)  (false, false)
     (false, false)  (false, false)  (false, false)



## Summary statistics

`summarize` function computes the following summary statistics of a `SnpArray`:  

* `maf`: minor allele frequencies, taking into account of missingness.  
* `minor_allele`: a `BitVector` indicating the minor allele for each SNP.   `minor_allele[j]==true` means A1 is the minor allele for SNP j; `minor_allele[j]==false` means A2 is the minor allele for SNP j.  
* `missings_by_snp`: number of missing genotypes for each snp.  
* `missings_by_person`: number of missing genotypes for each person.  


```julia
maf, minor_allele, missings_by_snp, missings_by_person = summarize(hapmap)
# minor allele frequencies
maf'
```




    1×13928 RowVector{Float64,Array{Float64,1}}:
     0.0  0.0776398  0.324074  0.191589  …  0.00154321  0.0417957  0.00617284




```julia
# total number of missing genotypes
sum(missings_by_snp), sum(missings_by_person)
```




    (11894, 11894)




```julia
# proportion of missing genotypes
sum(missings_by_snp) / length(hapmap)
```




    0.0026356890108565393



## Filtering

In almost all analyses, SNPs and individuals with low genotyping success rates are ignored. This filtering step is an important tool for removing likely false positives from association testing, as genotyping failure often occurs preferentially in cases or controls, or is correlated with the quantitative trait. `filter(s, min_success_rate_per_snp, min_success_rate_per_person)` does filtering according to the specified success rates for SNPs and people. Default is 0.98 for both.


```julia
# filtering SNPs and people to have both success rates above 0.98
snp_idx, person_idx = filter(hapmap, 0.98, 0.98)
# summary statistics of the filtered SnpArray
_, _, missings_by_snp_filtered, missings_by_person_filtered = summarize(hapmap[person_idx, snp_idx]);
```


```julia
# minimum SNP genotyping success rate after filtering ≥ 0.98
1.0 - maximum(missings_by_snp_filtered) / length(missings_by_person_filtered)
```




    0.9813084112149533




```julia
# minimum person genotyping success rate after filtering ≥ 0.98
1.0 - maximum(missings_by_person_filtered) / length(missings_by_snp_filtered)
```




    0.9818511796733213



## Random genotypes generation

`randgeno(a1freq)` generates a random genotype according to A1 allele frequency `a1freq`.


```julia
randgeno(0.5)
```




    (true, true)



`randgeno(maf, minor_allele)` generates a random genotype according to minor allele frequency `maf` and whether the minor allele is A1 (`minor_allele==true`) or A2 (`minor_allele==false`).


```julia
randgeno(0.25, true)
```




    (false, true)



`randgeno(n, maf, minor_allele)` generates a vector of random genotypes according to a common minor allele frequency `maf` and the minor allele.


```julia
randgeno(10, 0.25, true)
```




    10-element SnpArrays.SnpArray{1}:
     (false, true) 
     (false, true) 
     (true, true)  
     (true, true)  
     (true, true)  
     (true, true)  
     (true, true)  
     (true, true)  
     (false, false)
     (true, true)  



`randgeno(m, n, maf, minor_allele)` generates a random $m$-by-$n$ `SnpArray` according to a vector of minor allele frequencies `maf` and a minor allele indicator vector. The lengths of both vectors should be `n`.


```julia
# this is a random replicate of the hapmap data
randgeno(size(hapmap), maf, minor_allele)
```




    324×13928 SnpArrays.SnpArray{2}:
     (true, true)  (true, true)    …  (true, true)   (true, true) 
     (true, true)  (true, true)       (true, true)   (true, true) 
     (true, true)  (false, true)      (false, true)  (true, true) 
     (true, true)  (false, false)     (true, true)   (true, true) 
     (true, true)  (false, true)      (true, true)   (true, true) 
     (true, true)  (true, true)    …  (true, true)   (true, true) 
     (true, true)  (true, true)       (true, true)   (true, true) 
     (true, true)  (true, true)       (true, true)   (true, true) 
     (true, true)  (true, true)       (true, true)   (true, true) 
     (true, true)  (true, true)       (true, true)   (true, true) 
     (true, true)  (false, true)   …  (true, true)   (true, true) 
     (true, true)  (true, true)       (true, true)   (true, true) 
     (true, true)  (true, true)       (true, true)   (true, true) 
     ⋮                             ⋱                              
     (true, true)  (true, true)       (true, true)   (true, true) 
     (true, true)  (true, true)       (true, true)   (false, true)
     (true, true)  (true, true)       (true, true)   (true, true) 
     (true, true)  (false, true)   …  (true, true)   (true, true) 
     (true, true)  (true, true)       (false, true)  (true, true) 
     (true, true)  (false, true)      (true, true)   (true, true) 
     (true, true)  (true, true)       (true, true)   (true, true) 
     (true, true)  (true, true)       (true, true)   (true, true) 
     (true, true)  (false, true)   …  (true, true)   (true, true) 
     (true, true)  (true, true)       (true, true)   (true, true) 
     (true, true)  (true, true)       (false, true)  (true, true) 
     (true, true)  (true, true)       (true, true)   (true, true) 



## Subsetting

Subsetting a `SnpArray` works the same way as subsetting any other arrays.


```julia
# genotypes of the 1st person
hapmap[1, :]
```




    13928-element SnpArrays.SnpArray{1}:
     (true, true)  
     (true, true)  
     (false, false)
     (true, true)  
     (true, true)  
     (true, true)  
     (false, true) 
     (false, true) 
     (true, true)  
     (false, true) 
     (true, true)  
     (true, true)  
     (false, false)
     ⋮             
     (false, true) 
     (false, true) 
     (true, true)  
     (false, true) 
     (false, true) 
     (false, true) 
     (false, true) 
     (false, true) 
     (false, true) 
     (true, true)  
     (true, true)  
     (true, true)  




```julia
# genotypes of the 5th SNP
hapmap[:, 5]
```




    324-element SnpArrays.SnpArray{1}:
     (true, true)  
     (true, true)  
     (false, true) 
     (false, true) 
     (true, true)  
     (false, false)
     (false, false)
     (true, true)  
     (true, true)  
     (true, true)  
     (true, true)  
     (true, true)  
     (false, true) 
     ⋮             
     (false, false)
     (true, true)  
     (false, true) 
     (true, true)  
     (true, true)  
     (true, true)  
     (true, true)  
     (true, true)  
     (false, true) 
     (true, true)  
     (true, true)  
     (true, true)  




```julia
# subsetting both persons and SNPs
hapmap[1:5, 5:10]
```




    5×6 SnpArrays.SnpArray{2}:
     (true, true)   (true, true)  (false, true)  …  (true, true)   (false, true)
     (true, true)   (true, true)  (true, true)      (true, true)   (false, true)
     (false, true)  (true, true)  (true, true)      (false, true)  (true, true) 
     (false, true)  (true, true)  (true, true)      (true, true)   (false, true)
     (true, true)   (true, true)  (true, true)      (true, true)   (false, true)




```julia
# filter out rare SNPs with MAF < 0.05
hapmap[:, maf .≥ 0.05]
```




    324×12085 SnpArrays.SnpArray{2}:
     (true, true)   (false, false)  …  (false, true)  (false, true)
     (false, true)  (false, true)      (true, true)   (true, true) 
     (true, true)   (false, true)      (true, true)   (true, true) 
     (true, true)   (false, true)      (false, true)  (false, true)
     (true, true)   (false, true)      (true, true)   (true, true) 
     (false, true)  (true, true)    …  (false, true)  (false, true)
     (true, true)   (true, true)       (true, true)   (true, true) 
     (true, true)   (false, false)     (true, true)   (true, true) 
     (true, true)   (false, true)      (true, true)   (true, true) 
     (true, true)   (false, true)      (false, true)  (false, true)
     (true, true)   (false, true)   …  (true, true)   (true, true) 
     (true, true)   (true, true)       (false, true)  (false, true)
     (true, true)   (false, false)     (false, true)  (false, true)
     ⋮                              ⋱                              
     (true, true)   (false, true)      (false, true)  (false, true)
     (true, true)   (false, false)     (false, true)  (false, true)
     (true, true)   (false, false)     (true, true)   (true, true) 
     (true, true)   (false, false)  …  (true, true)   (true, true) 
     (true, true)   (false, true)      (true, true)   (true, true) 
     (true, true)   (true, true)       (false, true)  (false, true)
     (true, true)   (false, true)      (false, true)  (false, true)
     (true, true)   (false, true)      (true, true)   (true, true) 
     (true, true)   (false, false)  …  (false, true)  (false, true)
     (true, true)   (false, true)      (false, true)  (false, true)
     (true, true)   (false, false)     (false, true)  (false, true)
     (true, true)   (false, false)     (true, true)   (true, true) 




```julia
# filter out individuals with genotyping success rate < 0.90
hapmap[missings_by_person / people .< 0.1, :]
```




    220×13928 SnpArrays.SnpArray{2}:
     (true, true)  (true, true)   (false, false)  …  (true, true)   (true, true)
     (true, true)  (false, true)  (false, true)      (false, true)  (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)     (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)   …  (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)     (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (false, true)  (false, true)   …  (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)
     ⋮                                            ⋱                             
     (true, true)  (true, true)   (false, false)     (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)     (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)  …  (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)     (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)     (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)     (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)   …  (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)     (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)     (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)     (true, true)   (true, true)



`view` creates views of subarray without copying data and improve efficiency in many calculations.


```julia
mafcommon, = summarize(view(hapmap, :, maf .≥ 0.05))
mafcommon'
```




    1×12085 RowVector{Float64,Array{Float64,1}}:
     0.0776398  0.324074  0.191589  …  0.310937  0.23913  0.23913  0.23913



## Assignment

It is possible to assign specific genotypes to a `SnpArray` entry.


```julia
hapmap[1, 1]
```




    (true, true)




```julia
hapmap[1, 1] = (false, true)
hapmap[1, 1]
```




    (false, true)




```julia
hapmap[1, 1] = NaN
hapmap[1, 1]
```




    (true, false)




```julia
hapmap[1, 1] = 2
hapmap[1, 1]
```




    (true, true)



Subsetted assignment such as `hapmap[:, 1] = Nan` is also valid.

## Convert, copy and imputation

In most analyses we convert a whole `SnpArray` or slices of it to numeric arrays (matrix of **minor allele counts**) for statistical analysis. Keep in mind that the storage of resultant data can be up to 32 fold larger than that of the original `SnpArray`. Fortunately, rich collection of data types in `Julia` allow us choose one that fits into memory. Below are estimates of memory usage for some common data types with `n` persons and `p` SNPs. Here MAF denotes the **average** minor allele frequencies.

* `SnpArray`: $0.25np$ bytes  
* `Matrix{Int8}`: $np$ bytes  
* `Matrix{Float16}`: $2np$ bytes  
* `Matrix{Float32}`: $4np$ bytes  
* `Matrix{Float64}`: $8np$ bytes  
* `SparseMatrixCSC{Float64,Int64}`: $16 \cdot \text{NNZ} + 8(p+1) \approx 16np(2\text{MAF}(1-\text{MAF})+\text{MAF}^2) + 8(p+1) = 16np \cdot \text{MAF}(2-\text{MAF}) + 8(p+1)$ bytes. When the average MAF=0.25, this is about $7np$ bytes. When MAF=0.025, this is about $0.8np$ bypes, 10 fold smaller than the `Matrix{Float64}` type.  
* `SparseMatrixCSC{Int8,UInt32}`: $5 \cdot \text{NNZ} + 4(p+1) \approx 5np(2\text{MAF}(1-\text{MAF})+\text{MAF}^2) + 4(p+1) = 5np \cdot \text{MAF}(2-\text{MAF}) + 4(p+1)$ bytes. When the average MAF=0.25, this is about $2.2np$ bytes. When MAF=0.08, this is about $0.8np$ bypes, 10 fold smaller than `Matrix{Float64}` type.  
* Two `SparseMatrixCSC{Bool,Int64}`: $2np \cdot \text{MAF} \cdot 9 + 16(p+1) = 18 np \cdot \text{MAF} + 16(p+1)$ bytes. When the average MAF=0.25, this is about $4.5np$ bytes. When MAF=0.045, this is about $0.8np$ bytes, 10 fold smaller than `Matrix{Float64}` type.  

To be concrete, consider 2 typical data sets:  
* COPD (GWAS): $n = 6670$ individuals, $p = 630998$ SNPs, average MAF is 0.2454.
* GAW19 (sequencing study): $n = 959$ individuals, $p = 8348674$ SNPs, average MAF is 0.085.  

| Data Type | COPD | GAW19 |  
|---|---:|---:|  
| `SnpArray` | 1.05GB | 2GB |  
| `Matrix{Float64}` | 33.67GB | 64.05GB |  
| `SparseMatrixCSC{Float64,Int64}` | 29GB | 20.82GB |  
| `SparseMatrixCSC{Bool,Int64}` | 18.6GB | 12.386GB |  

Apparently for data sets with a majority of rare variants, converting to sparse matrices saves memory and often brings computational advantages too. In the `SparseMatrixCSC` format, the integer type of the row indices `rowval` and column pointer `colptr` should have maximal allowable value larger than the number of nonzeros in the matrix. The `InexactError()` error encountered during conversion often indicates that the integer type has a too small range. The utility function `estimatesize` conveniently estimates memory usage in bytes for the input data type.


```julia
# estimated memory usage if convert to Matrix{Float64}
estimatesize(people, snps, Matrix{Float64})
```




    3.6101376e7




```julia
# convert to Matrix{Float64}
hapmapf64 = convert(Matrix{Float64}, hapmap)
```




    324×13928 Array{Float64,2}:
     0.0  0.0  2.0  0.0  0.0  0.0  1.0  1.0  …  1.0  1.0  1.0  1.0  0.0  0.0  0.0
     0.0  1.0  1.0  1.0  0.0  0.0  0.0  1.0     0.0  0.0  0.0  0.0  0.0  1.0  0.0
     0.0  0.0  1.0  1.0  1.0  0.0  0.0  2.0     1.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.0  0.0  1.0  0.0  1.0  0.0  0.0  1.0     1.0  1.0  1.0  1.0  0.0  0.0  0.0
     0.0  0.0  1.0  1.0  0.0  0.0  0.0  2.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.0  1.0  0.0  0.0  2.0  0.0  0.0  0.0  …  2.0  1.0  1.0  1.0  0.0  0.0  0.0
     0.0  0.0  0.0  0.0  2.0  0.0  0.0  2.0     1.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.0  0.0  2.0  0.0  0.0  0.0  0.0  1.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.0  0.0  1.0  1.0  0.0  0.0  0.0  2.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.0  0.0  1.0  0.0  0.0  0.0  0.0  2.0     1.0  1.0  1.0  1.0  0.0  0.0  0.0
     0.0  0.0  1.0  0.0  0.0  0.0  0.0  2.0  …  0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.0  0.0  0.0  1.0  0.0  0.0  0.0  2.0     1.0  1.0  1.0  1.0  0.0  0.0  0.0
     0.0  0.0  2.0  0.0  1.0  0.0  0.0  1.0     2.0  1.0  1.0  1.0  0.0  0.0  0.0
     ⋮                        ⋮              ⋱                      ⋮            
     0.0  0.0  1.0  0.0  2.0  0.0  0.0  1.0     2.0  1.0  1.0  1.0  0.0  0.0  0.0
     0.0  0.0  2.0  0.0  0.0  0.0  0.0  2.0     2.0  1.0  1.0  1.0  0.0  0.0  0.0
     0.0  0.0  2.0  0.0  1.0  0.0  0.0  1.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.0  0.0  2.0  1.0  0.0  0.0  0.0  1.0  …  0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.0  0.0  1.0  0.0  0.0  0.0  0.0  1.0     1.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0     2.0  1.0  1.0  1.0  0.0  0.0  0.0
     0.0  0.0  1.0  1.0  0.0  0.0  0.0  2.0     1.0  1.0  1.0  1.0  0.0  0.0  0.0
     0.0  0.0  1.0  0.0  0.0  0.0  0.0  2.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.0  0.0  2.0  0.0  1.0  0.0  0.0  1.0  …  1.0  1.0  1.0  1.0  0.0  0.0  0.0
     0.0  0.0  1.0  1.0  0.0  0.0  0.0  1.0     1.0  1.0  1.0  1.0  0.0  0.0  0.0
     0.0  0.0  2.0  1.0  0.0  0.0  0.0  2.0     1.0  1.0  1.0  1.0  0.0  0.0  0.0
     0.0  0.0  2.0  0.0  0.0  0.0  0.0  2.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0




```julia
# actual memory usage of Matrix{Float64}
Base.summarysize(hapmapf64)
```




    36101376




```julia
# average maf of the hapmap3 data set
mean(maf)
```




    0.222585591341583




```julia
# estimated memory usage if convert to SparseMatrixCSC{Float32, UInt32} matrix
estimatesize(people, snps, SparseMatrixCSC{Float32, UInt32}, mean(maf))
```




    1.4338389205819245e7




```julia
# convert to SparseMatrixCSC{Float32, UInt32} matrix
hapmapf32sp = convert(SparseMatrixCSC{Float32, UInt32}, hapmap)
```




    324×13928 SparseMatrixCSC{Float32,UInt32} with 1614876 stored entries:
      [2    ,     2]  =  1.0
      [6    ,     2]  =  1.0
      [15   ,     2]  =  1.0
      [31   ,     2]  =  1.0
      [33   ,     2]  =  1.0
      [35   ,     2]  =  1.0
      [43   ,     2]  =  1.0
      [44   ,     2]  =  1.0
      [50   ,     2]  =  1.0
      [54   ,     2]  =  1.0
      ⋮
      [135  , 13927]  =  1.0
      [148  , 13927]  =  1.0
      [160  , 13927]  =  1.0
      [164  , 13927]  =  2.0
      [167  , 13927]  =  1.0
      [185  , 13927]  =  1.0
      [266  , 13927]  =  1.0
      [280  , 13927]  =  1.0
      [288  , 13927]  =  1.0
      [118  , 13928]  =  2.0
      [231  , 13928]  =  2.0




```julia
# actual memory usage if convert to SparseMatrixCSC{Float32, UInt32} matrix
Base.summarysize(hapmapf32sp)
```




    12974764



By default the `convert()` method converts missing genotypes to `NaN`.


```julia
# number of missing genotypes
countnz(isnan.(hapmap)), countnz(isnan.(hapmapf64))
```




    (11894, 11894)



One can enforce **crude imputation** by setting the optional argument `impute=true`. Imputation is done by generating two random alleles according to the minor allele frequency. This is a neutral but not an optimal strategy, and users should impute missing genotypes by more advanced methods.


```julia
hapmapf64impute = convert(Matrix{Float64}, hapmap; impute = true)
countnz(isnan.(hapmapf64impute))
```




    0



By default `convert()` translates genotypes according to the *additive* SNP model, which essentially counts the number of **minor allele** (0, 1 or 2) per genotype. Other SNP models are *dominant* and *recessive*, both in terms of the **minor allele**. When `A1` is the minor allele, genotypes are translated to real number according to

| Genotype | `SnpArray` | `model=:additive` | `model=:dominant` | `model=:recessive` |    
|:---:|:---:|:---:|:---:|:---:|  
| A1,A1 | (false,false) | 2 | 1 | 1 |  
| A1,A2 | (false,true) | 1 | 1 | 0 |  
| A2,A2 | (true,true) | 0 | 0 | 0 |  
| missing | (true,false) | NaN | NaN | NaN | 

When `A2` is the minor allele, genotypes are translated according to

| Genotype | `SnpArray` | `model=:additive` | `model=:dominant` | `model=:recessive` |    
|:---:|:---:|:---:|:---:|:---:|  
| A1,A1 | (false,false) | 0 | 0 | 0 |  
| A1,A2 | (false,true) | 1 | 1 | 0 |  
| A2,A2 | (true,true) | 2 | 1 | 1 |  
| missing | (true,false) | NaN | NaN | NaN |


```julia
[convert(Vector{Float64}, hapmap[1:10, 5]; model = :additive) convert(Vector{Float64}, hapmap[1:10, 5]; model = :dominant) convert(Vector{Float64}, hapmap[1:10, 5]; model = :recessive)]
```




    10×3 Array{Float64,2}:
     0.0  0.0  0.0
     0.0  0.0  0.0
     1.0  1.0  0.0
     1.0  1.0  0.0
     0.0  0.0  0.0
     2.0  1.0  1.0
     2.0  1.0  1.0
     0.0  0.0  0.0
     0.0  0.0  0.0
     0.0  0.0  0.0



By default `convert()` does **not** center and scale genotypes. Setting the optional arguments `center=true, scale=true` centers genotypes at 2MAF and scales them by $[2 \cdot \text{MAF} \cdot (1 - \text{MAF})]^{-1/2}$. Mono-allelic SNPs (MAF=0) are not scaled.


```julia
[convert(Vector{Float64}, hapmap[:, 5]) convert(Vector{Float64}, hapmap[:, 5]; center = true, scale = true)]
```




    324×2 Array{Float64,2}:
     0.0  -1.25702 
     0.0  -1.25702 
     1.0   0.167017
     1.0   0.167017
     0.0  -1.25702 
     2.0   1.59106 
     2.0   1.59106 
     0.0  -1.25702 
     0.0  -1.25702 
     0.0  -1.25702 
     0.0  -1.25702 
     0.0  -1.25702 
     1.0   0.167017
     ⋮             
     2.0   1.59106 
     0.0  -1.25702 
     1.0   0.167017
     0.0  -1.25702 
     0.0  -1.25702 
     0.0  -1.25702 
     0.0  -1.25702 
     0.0  -1.25702 
     1.0   0.167017
     0.0  -1.25702 
     0.0  -1.25702 
     0.0  -1.25702 



`copy!()` is the in-place version of `convert()`. Options such as GWAS loop over SNPs and perform statistical anlaysis for each SNP. This can be achieved by


```julia
g = zeros(people)
for j = 1:snps
    copy!(g, hapmap[:, j]; model = :additive, impute = true)
    # do statistical anlaysis
end
```

## Empirical kinship matrix

`grm` function computes the empirical kinship matrix using either the genetic relationship matrix, `grm(A, model=:GRM)`, or the method of moment method, `grm(A, model=:MoM)`. 

!!! note

    Missing genotypes are imputed according to minor allele frequencies on the fly.  
    


!!! note

    It is often necessary to filter SNPs according to minor allele frequency and LD before calculating empirical kinship matrix.  


By default, `grm` exlcude SNPs with minor allele frequency below 0.01. This can be changed by the keyword argument `maf_threshold`.


```julia
# GRM using all SNPs with MAF ≥ 0.01. Same as
# grm(hapmap; maf_threshold = 0.01)
grm(hapmap)
```




    324×324 Array{Float64,2}:
     0.571688   0.0462807  0.0192437  …  0.0650573  0.0716031  0.0653311
     0.0462807  0.544895   0.0286567     0.051311   0.0449147  0.0629484
     0.0192437  0.0286567  0.520718      0.0460941  0.0297353  0.0362593
     0.0486305  0.0369559  0.0293866     0.0599803  0.0659248  0.0599897
     0.052624   0.0423282  0.0247818     0.0715649  0.0581852  0.0658928
     0.0441266  0.0319019  0.0378722  …  0.0697301  0.0563314  0.065507 
     0.0397138  0.0217914  0.0120808     0.0446991  0.0381401  0.0371412
     0.0410787  0.0382596  0.021239      0.0577011  0.0547265  0.065843 
     0.030168   0.0315418  0.0169194     0.0334515  0.0460561  0.0380265
     0.0394529  0.0422534  0.0258607     0.0665238  0.0574226  0.0496232
     0.0477793  0.0454803  0.0230092  …  0.0585435  0.0630797  0.0610807
     0.0607075  0.0387026  0.0370969     0.0669321  0.0590524  0.0709323
     0.0357632  0.0438016  0.0264945     0.0571821  0.0701542  0.0635772
     ⋮                                ⋱                                 
     0.0583693  0.0556652  0.0382274     0.0864722  0.0864521  0.0971764
     0.0676922  0.0521073  0.0324018     0.091956   0.0803552  0.0859155
     0.0679075  0.0558743  0.0386209     0.0899731  0.0844895  0.0857414
     0.0652022  0.0564352  0.0390764  …  0.0878313  0.0765041  0.0892029
     0.063271   0.0561774  0.0364754     0.0789625  0.0734351  0.0822378
     0.0648145  0.0612811  0.040045      0.0896104  0.0901562  0.0812354
     0.0625796  0.0600032  0.0380649     0.0831302  0.0852875  0.0809756
     0.0622714  0.0629473  0.0334698     0.0932245  0.0900163  0.0869392
     0.0660527  0.0633805  0.0322493  …  0.0973118  0.0791239  0.0866666
     0.0650573  0.051311   0.0460941     0.60326    0.0819594  0.0901512
     0.0716031  0.0449147  0.0297353     0.0819594  0.595108   0.0823624
     0.0653311  0.0629484  0.0362593     0.0901512  0.0823624  0.582986 




```julia
# GRM using all SNPs with MAF ≥ 0.05
grm(hapmap; maf_threshold = 0.05)
```




    324×324 Array{Float64,2}:
     0.569625   0.0535704  0.0232651  …  0.0737173  0.077806   0.0728392
     0.0535704  0.542874   0.0292732     0.0564599  0.0502244  0.0685582
     0.0232651  0.0292732  0.52172       0.0518046  0.0350828  0.0406411
     0.0530896  0.0398713  0.0315497     0.0621363  0.0726467  0.0656261
     0.0552446  0.0461061  0.0260241     0.0764474  0.0638611  0.0722666
     0.0491941  0.0332782  0.0418767  …  0.075888   0.0612478  0.0647859
     0.0426535  0.0262393  0.0149532     0.0467589  0.0381211  0.043397 
     0.044185   0.0420889  0.0242141     0.0573316  0.0569442  0.0666711
     0.0312286  0.0339761  0.0202836     0.0375401  0.0499576  0.041732 
     0.0418681  0.0447897  0.0290518     0.0711012  0.0638013  0.0531209
     0.0551049  0.047964   0.0237454  …  0.0628456  0.068017   0.0672129
     0.0663155  0.042911   0.0417686     0.0743179  0.0654391  0.0735138
     0.0359789  0.0473879  0.0296654     0.0629452  0.0732141  0.0676924
     ⋮                                ⋱                                 
     0.0609321  0.05769    0.0424677     0.0922023  0.0921072  0.105169 
     0.0720513  0.0550647  0.0373036     0.0979529  0.0860333  0.0926531
     0.0673758  0.0596123  0.0407195     0.0922654  0.090949   0.0909565
     0.0706774  0.064258   0.0456276  …  0.0946531  0.0831749  0.093832 
     0.0700268  0.0596294  0.0393517     0.0864101  0.0816786  0.0898796
     0.0692969  0.0634394  0.0440346     0.0938751  0.0951635  0.0888866
     0.0663192  0.0625808  0.0413768     0.0858972  0.0889115  0.0878044
     0.068348   0.067889   0.039336      0.100281   0.0920632  0.0865986
     0.0721291  0.0674194  0.0388886  …  0.106512   0.0845076  0.0901555
     0.0737173  0.0564599  0.0518046     0.59046    0.0903191  0.0945729
     0.077806   0.0502244  0.0350828     0.0903191  0.588007   0.0883015
     0.0728392  0.0685582  0.0406411     0.0945729  0.0883015  0.570052 




```julia
# GRM using every other SNP, with maf ≥ 0.01
grm(view(hapmap, :, 1:2:snps))
```




    324×324 Array{Float64,2}:
     0.55891    0.042845   0.027682   …  0.0677721  0.0746666  0.0687635
     0.042845   0.558745   0.0286764     0.0578987  0.0462962  0.0569583
     0.027682   0.0286764  0.512307      0.0390934  0.0396354  0.0476762
     0.045369   0.0469047  0.028284      0.0515663  0.0621996  0.0564942
     0.0523796  0.0477215  0.0267944     0.0676464  0.0566912  0.0656048
     0.0523564  0.03959    0.040644   …  0.0775905  0.0623572  0.0523832
     0.0393927  0.0262359  0.0178367     0.0471618  0.0347238  0.0343188
     0.047151   0.0379639  0.0263928     0.0590008  0.0559667  0.0688052
     0.0269253  0.0227214  0.0201293     0.0331204  0.0438332  0.0359987
     0.0323098  0.0402429  0.0238774     0.0598208  0.0408401  0.0505469
     0.0494403  0.0504749  0.0211352  …  0.0633784  0.0648717  0.0534096
     0.0626913  0.0494432  0.0407554     0.0668062  0.0648524  0.0606319
     0.0256349  0.0433059  0.0229496     0.0547733  0.0687123  0.063634 
     ⋮                                ⋱                                 
     0.0605549  0.0506841  0.0403012     0.0871397  0.0847116  0.101929 
     0.0663559  0.0599439  0.0310283     0.0894028  0.0762304  0.0813586
     0.0694759  0.0543091  0.0394317     0.0846661  0.0820957  0.096536 
     0.0654772  0.0581756  0.0404721  …  0.0934607  0.0768118  0.0828469
     0.0646551  0.0613316  0.0423482     0.0735479  0.0703157  0.0875644
     0.066537   0.0632911  0.0413394     0.0763028  0.082904   0.0844916
     0.0638155  0.0649174  0.032502      0.0869326  0.0919846  0.0835989
     0.0598241  0.0613027  0.0346562     0.0946688  0.0910173  0.0840148
     0.076463   0.0668378  0.0405983  …  0.101928   0.0747384  0.0896328
     0.0677721  0.0578987  0.0390934     0.573711   0.0700193  0.0959231
     0.0746666  0.0462962  0.0396354     0.0700193  0.57723    0.0752888
     0.0687635  0.0569583  0.0476762     0.0959231  0.0752888  0.588204 




```julia
# MoM using all SNPs with MAF ≥ 0.01
grm(hapmap; method = :MoM)
```




    324×324 Array{Float64,2}:
     0.539321    0.0346778  0.00253189  …  0.0537053  0.0637509  0.0503962
     0.0346778   0.517693   0.0141139      0.0421233  0.0395233  0.0499234
     0.00253189  0.0141139  0.499493       0.0320777  0.0210867  0.0188412
     0.0436597   0.0291232  0.0233322      0.0526417  0.0696601  0.0496871
     0.0456688   0.0333778  0.0156502      0.0658782  0.0570144  0.0564235
     0.0317232   0.0197867  0.0261686   …  0.0587872  0.0496871  0.0488598
     0.0243958   0.0104502  0.00276825     0.0313686  0.0305414  0.0238049
     0.0256958   0.0285322  0.00879561     0.0384596  0.0460234  0.0462598
     0.0220322   0.0260504  0.0135229      0.0298323  0.0488598  0.0301868
     0.0209685   0.0255777  0.0104502      0.0526417  0.0508689  0.0319596
     0.0362142   0.0330232  0.00525011  …  0.0496871  0.0561872  0.0500416
     0.0483871   0.0313686  0.0268777      0.0546508  0.0546508  0.051578 
     0.0313686   0.0415324  0.0238049      0.0585508  0.0704874  0.0577235
     ⋮                                  ⋱                                 
     0.0416506   0.041887   0.0254595      0.0713146  0.0792329  0.081833 
     0.0620963   0.0482689  0.0267595      0.0917604  0.0792329  0.0836057
     0.0552417   0.0453143  0.0256958      0.0806511  0.0817148  0.0727328
     0.0454325   0.0421233  0.0235685   …  0.072851   0.0628054  0.0703692
     0.0558326   0.0534689  0.0297141      0.0775784  0.070251   0.0736783
     0.0520507   0.0542962  0.0312505      0.0788784  0.0866785  0.0714328
     0.0522871   0.0507507  0.0295959      0.0772238  0.0827784  0.0732056
     0.0482689   0.0521689  0.0196685      0.0820693  0.0838421  0.0720237
     0.0500416   0.0477961  0.0214412   …  0.0880967  0.0735601  0.0693055
     0.0537053   0.0421233  0.0320777      0.561303   0.0807693  0.0755692
     0.0637509   0.0395233  0.0210867      0.0807693  0.566858   0.0688328
     0.0503962   0.0499234  0.0188412      0.0755692  0.0688328  0.533766 



## Principal components 

Principal compoenent analysis is widely used in genome-wide association analysis (GWAS) for adjusting population substructure. `pca(A, pcs)` computes the top `pcs` principal components of a `SnpArray`. Each SNP is centered at $2\text{MAF}$ and scaled by $[2\text{MAF}(1-\text{MAF})]^{-1/2}$. The output is  

* `pcscore`: top `pcs` eigen-SNPs, or principal scores, in each column  
* `pcloading`: top `pcs` eigen-vectors, or principal loadings, in each column  
* `pcvariance`: top `pcs` eigen-values, or principal variances

Missing genotypes are imputed according the minor allele frequencies on the fly. This implies that, in the presence of missing genotypes, running the function on the same `SnpArray` twice may produce slightly different answers. For reproducibility, it is a good practice to set the random seed before each function that does imputation on the fly.


```julia
srand(123) # set seed
pcscore, pcloading, pcvariance = pca(hapmap, 3)
```




    ([-38.7231 -1.2983 -7.00541; -32.6096 -1.21052 -3.3232; … ; -48.9263 -2.06102 2.17374; -48.8627 0.274894 6.49518], [2.56616e-19 8.19569e-19 5.52006e-19; 0.00143962 -0.0042375 -0.00311816; … ; 0.00313326 -0.00427486 -0.0152038; -9.09523e-5 -0.00287777 0.0037855], [1841.4, 225.324, 70.7084])



To use eigen-SNPs for plotting or as covariates in GWAS, we typically scale them by their standard deviations so that they have mean zero and unit variance.


```julia
# standardize eigen-SNPs before plotting or GWAS
scale!(pcscore, 1 ./ sqrt.(pcvariance))
std(pcscore, 1)
```




    1×3 Array{Float64,2}:
     1.0  1.0  1.0



Internally `pca` converts `SnpArray` to the matrix of minor allele counts. The default format is `Matrix{Float64}`, which can easily exceed memory limit. Users have several options when the default `Matrix{Float64}` cannot fit into memory.  

* Use other intermediate matrix types.


```julia
# use single precision matrix and display the principal variances
# approximately same answer as double precision
srand(123)
pca(hapmap, 3, Matrix{Float32})[3]
```




    3-element Array{Float32,1}:
     1841.39  
      225.324 
       70.7084



* Use subset of SNPs


```julia
# principal components using every other SNP capture about half the variance
srand(123)
pca(view(hapmap, :, 1:2:snps), 3)[3]
```




    3-element Array{Float64,1}:
     926.622 
     113.188 
      36.4866



* Use sparse matrix. For large data sets with majority of rare variants, `pca_sp` is more efficient by first converting `SnpArray` to a sparse matrix (default is `SparseMatrixCSC{Float64, Int64}`) and then computing principal components using iterative algorithms. 


```julia
# approximately same answer if we use Float16 sparse matrix
srand(123)
pca_sp(hapmap, 3, SparseMatrixCSC{Float32, UInt32})[3]
```




    3-element Array{Float32,1}:
     1841.39  
      225.324 
       70.7084


