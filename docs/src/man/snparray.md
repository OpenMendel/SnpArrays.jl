
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



* `SnpArray` can be initialized from a matrix of A2 allele counts.


```julia
SnpArray(rand(0:2, 5, 3))
```




    5×3 SnpArrays.SnpArray{2}:
     (false, true)   (true, true)    (false, true)
     (false, true)   (true, true)    (true, true) 
     (false, true)   (true, true)    (true, true) 
     (false, false)  (false, false)  (true, true) 
     (true, true)    (false, true)   (true, true) 



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
maf
```




    13928-element Array{Float64,1}:
     0.0       
     0.0776398 
     0.324074  
     0.191589  
     0.441358  
     0.0       
     0.00462963
     0.453704  
     0.226852  
     0.14486   
     0.0       
     0.483025  
     0.25387   
     ⋮         
     0.239938  
     0.239938  
     0.255486  
     0.23913   
     0.238318  
     0.310937  
     0.23913   
     0.23913   
     0.23913   
     0.00154321
     0.0417957 
     0.00617284




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
1 - maximum(missings_by_snp_filtered) / length(missings_by_person_filtered)
```




    0.9813084112149533




```julia
# minimum individual genotyping success rate after filtering ≥ 0.98
1 - maximum(missings_by_person_filtered) / length(missings_by_snp_filtered)
```




    0.9818511796733213



## Random genotypes generation

`randgeno(a1freq)` generates a random genotype according to A1 allele frequency `a1freq`.


```julia
randgeno(0.5)
```




    (false, false)



`randgeno(maf, minor_allele)` generates a random genotype according to minor allele frequency `maf` and whether the minor allele is A1 (`minor_allele==true`) or A2 (`minor_allele==false`).


```julia
randgeno(0.25, true)
```




    (true, true)



`randgeno(n, maf, minor_allele)` generates a vector of random genotypes according to a common minor allele frequency `maf` and the minor allele.


```julia
randgeno(10, 0.25, true)
```




    10-element SnpArrays.SnpArray{1}:
     (false, true) 
     (false, true) 
     (false, true) 
     (false, true) 
     (true, true)  
     (true, true)  
     (false, false)
     (true, true)  
     (false, true) 
     (false, true) 



`randgeno(m, n, maf, minor_allele)` generates a random $m$-by-$n$ `SnpArray` according to a vector of minor allele frequencies `maf` and a minor allele indicator vector. The lengths of both vectors should be `n`.


```julia
# this is a random replicate of the hapmap data
randgeno(size(hapmap), maf, minor_allele)
```




    324×13928 SnpArrays.SnpArray{2}:
     (true, true)  (true, true)   (true, true)   …  (true, true)   (true, true)
     (true, true)  (false, true)  (false, true)     (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)      (false, true)  (true, true)
     (true, true)  (true, true)   (true, true)      (true, true)   (true, true)
     (true, true)  (false, true)  (true, true)      (true, true)   (true, true)
     (true, true)  (false, true)  (false, true)  …  (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)     (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)     (true, true)   (true, true)
     (true, true)  (false, true)  (false, true)     (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)  …  (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)      (true, true)   (true, true)
     ⋮                                           ⋱                             
     (true, true)  (false, true)  (false, true)     (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)     (true, true)   (true, true)
     (true, true)  (false, true)  (true, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)   …  (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)     (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)     (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)     (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)   …  (true, true)   (true, true)
     (true, true)  (false, true)  (true, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)     (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)      (true, true)   (true, true)



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
mafcommon
```




    12085-element Array{Float64,1}:
     0.0776398
     0.324074 
     0.191589 
     0.441358 
     0.453704 
     0.226852 
     0.14486  
     0.483025 
     0.25387  
     0.109907 
     0.221875 
     0.475232 
     0.305556 
     ⋮        
     0.253125 
     0.238318 
     0.235016 
     0.239938 
     0.239938 
     0.255486 
     0.23913  
     0.238318 
     0.310937 
     0.23913  
     0.23913  
     0.23913  



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
     0.571458   0.0464311  0.0199551  …  0.065172   0.0719459  0.0653116
     0.0464311  0.544725   0.0282826     0.0512752  0.0449652  0.0631871
     0.0199551  0.0282826  0.520116      0.046552   0.0305444  0.0366185
     0.0485422  0.0370034  0.0294222     0.0600245  0.0658291  0.0599777
     0.0527247  0.0428915  0.0250871     0.0716037  0.0587246  0.0657509
     0.0449014  0.0324916  0.0392128  …  0.0705846  0.0570837  0.0653063
     0.0395131  0.0224096  0.0118297     0.0444225  0.0383403  0.0372241
     0.0414288  0.0385228  0.0216014     0.0578062  0.0549297  0.0663108
     0.0302047  0.0311982  0.0167017     0.03367    0.04633    0.0380357
     0.0391444  0.0420736  0.0262111     0.0667603  0.057535   0.048993 
     0.0483412  0.0453938  0.0229784  …  0.0586796  0.0629187  0.0611821
     0.0607361  0.0386021  0.0373459     0.0661212  0.0586743  0.0698014
     0.036619   0.0441432  0.0266479     0.0574166  0.06993    0.0636137
     ⋮                                ⋱                                 
     0.057858   0.0550421  0.0385672     0.0862127  0.0867892  0.0968649
     0.0670946  0.051947   0.0329633     0.0918435  0.0803133  0.0861184
     0.06812    0.056159   0.038998      0.0898582  0.0848113  0.0854561
     0.0654276  0.0568588  0.0393993  …  0.0874573  0.0767071  0.089505 
     0.063116   0.056254   0.0357547     0.0785316  0.073411   0.0817352
     0.064895   0.0609217  0.0399485     0.0886267  0.0902875  0.0813491
     0.0623362  0.0601209  0.0385246     0.0822759  0.0841088  0.0810898
     0.0618955  0.0630099  0.0332965     0.0933279  0.0898342  0.0863454
     0.066011   0.063046   0.0324327  …  0.0972288  0.0795694  0.0866725
     0.065172   0.0512752  0.046552      0.603114   0.0821105  0.0903747
     0.0719459  0.0449652  0.0305444     0.0821105  0.595455   0.0825498
     0.0653116  0.0631871  0.0366185     0.0903747  0.0825498  0.582903 




```julia
# GRM using all SNPs with MAF ≥ 0.05
grm(hapmap; maf_threshold = 0.05)
```




    324×324 Array{Float64,2}:
     0.569655   0.0533259  0.0238227  …  0.0733238  0.0773046  0.0725078
     0.0533259  0.543035   0.0292765     0.0565683  0.0502522  0.0689859
     0.0238227  0.0292765  0.521168      0.0523316  0.0353608  0.0407409
     0.0530866  0.0402229  0.0320051     0.0619395  0.072497   0.0651918
     0.0558751  0.045637   0.0260024     0.0764646  0.0639329  0.0724448
     0.0487862  0.0324028  0.0414961  …  0.075731   0.0607731  0.0646571
     0.0425807  0.0259461  0.0158766     0.0469247  0.0383466  0.0437265
     0.0445209  0.0415715  0.0246058     0.0574699  0.0568312  0.0665313
     0.0313927  0.0335657  0.01984       0.0369867  0.0502879  0.0415424
     0.0419546  0.0443697  0.02947       0.0709567  0.0633844  0.0528883
     0.0547824  0.0481855  0.024218   …  0.0626619  0.0685641  0.0672536
     0.0668034  0.0429604  0.0412379     0.0736168  0.0656362  0.0739072
     0.0361395  0.0472199  0.0295056     0.0632602  0.0731142  0.0671983
     ⋮                                ⋱                                 
     0.0604711  0.0578328  0.0430631     0.0922183  0.092174   0.105441 
     0.0719807  0.0549142  0.0378689     0.0980107  0.0858935  0.0927524
     0.0675168  0.0597127  0.0417714     0.0923096  0.0912273  0.0915402
     0.0705962  0.0636391  0.0453237  …  0.0942765  0.0835207  0.093857 
     0.0698496  0.0590258  0.0398921     0.0863023  0.081476   0.0898919
     0.0695019  0.0642225  0.0446517     0.0941713  0.0955554  0.088804 
     0.0672219  0.0627248  0.040251      0.0852741  0.0894345  0.0869601
     0.069019   0.0682414  0.0396043     0.0996116  0.092133   0.0865049
     0.0722643  0.0673632  0.0387407  …  0.106599   0.0842294  0.0901445
     0.0733238  0.0565683  0.0523316     0.590262   0.0897775  0.0944278
     0.0773046  0.0502522  0.0353608     0.0897775  0.587805   0.0880512
     0.0725078  0.0689859  0.0407409     0.0944278  0.0880512  0.570448 




```julia
# GRM using every other SNP, with maf ≥ 0.01
grm(view(hapmap, :, 1:2:snps))
```




    324×324 Array{Float64,2}:
     0.558675   0.0431596  0.0281906  …  0.0678361  0.0737569  0.0678428
     0.0431596  0.557339   0.0288423     0.0584065  0.0460293  0.0559139
     0.0281906  0.0288423  0.512262      0.0404619  0.0395088  0.047475 
     0.0454521  0.0464444  0.0267929     0.0518465  0.0617595  0.0566491
     0.0523056  0.0499494  0.0265126     0.0678693  0.0567765  0.0654963
     0.0524025  0.0409955  0.0405314  …  0.0769924  0.0626306  0.0530175
     0.0396371  0.0266048  0.0170996     0.0472069  0.0345508  0.0343661
     0.0473392  0.0388329  0.0270152     0.059386   0.0558696  0.0688371
     0.0266013  0.0237236  0.0193008     0.033394   0.043639   0.0361471
     0.0320554  0.0396779  0.0241781     0.0599865  0.0410583  0.0504274
     0.0491539  0.0504333  0.0210508  …  0.0629897  0.0652705  0.0533614
     0.0620909  0.0497012  0.0414584     0.0667277  0.0650775  0.0609132
     0.0250203  0.0432866  0.0219478     0.0547627  0.0687206  0.0634917
     ⋮                                ⋱                                 
     0.0610989  0.0504615  0.0411155     0.0871821  0.0857404  0.102135 
     0.0667531  0.0606607  0.03117       0.0898112  0.0760716  0.0812821
     0.0695847  0.0541784  0.0394035     0.0848243  0.0817369  0.0958491
     0.0656943  0.0586907  0.0406232  …  0.0935051  0.0763204  0.0827709
     0.0639221  0.0619239  0.0409911     0.0730449  0.069307   0.087749 
     0.0666386  0.0638029  0.0420165     0.0764784  0.0830993  0.0850801
     0.0639109  0.0642512  0.0320533     0.0867362  0.0932944  0.0844009
     0.0604671  0.0619833  0.0347137     0.09504    0.0918656  0.0838407
     0.0766471  0.0665513  0.0397895  …  0.102109   0.0746396  0.0900517
     0.0678361  0.0584065  0.0404619     0.573673   0.0698129  0.0958643
     0.0737569  0.0460293  0.0395088     0.0698129  0.576771   0.0746628
     0.0678428  0.0559139  0.047475      0.0958643  0.0746628  0.588228 




```julia
# MoM using all SNPs with MAF ≥ 0.01
grm(hapmap; method = :MoM)
```




    324×324 Array{Float64,2}:
     0.539321    0.0355051  0.00253189  …  0.0538235  0.0633964  0.0511053
     0.0355051   0.518284   0.0147048      0.0423597  0.0397597  0.0508689
     0.00253189  0.0147048  0.499257       0.0324323  0.0214412  0.0201412
     0.0434234   0.0292413  0.0235685      0.0526417  0.0694237  0.0500416
     0.0454325   0.0339687  0.0162412      0.0662328  0.057369   0.0568963
     0.0317232   0.0204958  0.0262868   …  0.0592599  0.0500416  0.0492143
     0.0249867   0.0113956  0.00359554     0.0317232  0.0310141  0.0252231
     0.0261686   0.0285322  0.00974108     0.0384596  0.0463779  0.0469689
     0.0217958   0.0266413  0.014823       0.0299504  0.0493325  0.0307777
     0.0215594   0.0256958  0.0109229      0.0529962  0.0508689  0.0325505
     0.0357414   0.0333778  0.00584102  …  0.0495689  0.056069   0.050278 
     0.0473234   0.0321959  0.027114       0.0550053  0.0550053  0.0518144
     0.031605    0.0425961  0.0254595      0.0587872  0.0706055  0.0580781
     ⋮                                  ⋱                                 
     0.041887    0.0420052  0.0264049      0.0713146  0.0789966  0.0824239
     0.0622145   0.0480325  0.027705       0.0917604  0.0792329  0.0836057
     0.054769    0.0459052  0.0265231      0.0807693  0.0815966  0.0727328
     0.0449597   0.0429506  0.0241595   …  0.0727328  0.0626873  0.070251 
     0.0564235   0.052878   0.0291232      0.0772238  0.0707237  0.0742692
     0.0533508   0.0551235  0.0313686      0.0794693  0.0865603  0.0719056
     0.0519326   0.050278   0.0285322      0.0771056  0.0827784  0.0733238
     0.0488598   0.0527598  0.0209685      0.0823057  0.0841966  0.0723783
     0.0506325   0.047678   0.0217958   …  0.087624   0.0734419  0.068951 
     0.0538235   0.0423597  0.0324323      0.561185   0.0811239  0.0763965
     0.0633964   0.0397597  0.0214412      0.0811239  0.566858   0.0691873
     0.0511053   0.0508689  0.0201412      0.0763965  0.0691873  0.533648 



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



## Linear algebra with SnpArray

In some applications we want to perform linear algebra using SnpArray directly without expanding it to numeric matrix. The implementation assumes:  
1. the SnpArray does not have missing genotypes, and  
2. the matrix corresponding to SnpArray is the matrix of A2 allele counts.


```julia
# generate a random SnpArray according to A2 allele counts
n, p = 100, 200
S = SnpArray(rand(0:2, n, p))
```




    100×200 SnpArrays.SnpArray{2}:
     (false, false)  (false, false)  …  (false, false)  (false, false)
     (true, true)    (false, true)      (false, false)  (false, false)
     (false, true)   (false, true)      (true, true)    (false, true) 
     (false, false)  (true, true)       (false, false)  (false, false)
     (true, true)    (false, true)      (true, true)    (false, true) 
     (false, false)  (true, true)    …  (false, false)  (false, false)
     (true, true)    (false, true)      (false, false)  (false, true) 
     (false, false)  (false, true)      (false, true)   (true, true)  
     (false, true)   (false, true)      (false, false)  (true, true)  
     (false, false)  (false, true)      (false, true)   (false, false)
     (false, true)   (false, true)   …  (true, true)    (false, false)
     (false, true)   (true, true)       (false, false)  (false, true) 
     (false, false)  (false, false)     (true, true)    (false, true) 
     ⋮                               ⋱                                
     (false, true)   (false, false)     (true, true)    (false, false)
     (true, true)    (false, false)     (false, true)   (true, true)  
     (false, false)  (false, true)   …  (true, true)    (false, false)
     (false, true)   (false, true)      (false, true)   (true, true)  
     (false, false)  (false, false)     (false, true)   (false, true) 
     (false, true)   (false, true)      (false, true)   (false, true) 
     (false, false)  (true, true)       (true, true)    (false, false)
     (true, true)    (true, true)    …  (false, true)   (false, true) 
     (false, true)   (false, false)     (false, true)   (true, true)  
     (false, true)   (false, false)     (false, false)  (false, true) 
     (true, true)    (false, false)     (false, false)  (true, true)  
     (false, false)  (false, true)      (false, true)   (false, false)



* (snparray)-vector multiplications


```julia
S * randn(p)
```




    100-element Array{Float64,1}:
     30.615   
     -7.4687  
     23.0195  
      0.817897
     11.5828  
     18.367   
      2.00545 
     23.3611  
     17.7711  
     15.1654  
     12.0407  
     16.3625  
     10.5148  
      ⋮       
     18.9264  
     18.157   
     31.7193  
     23.6554  
     15.2196  
     25.8444  
     12.3667  
     12.3029  
     22.1251  
      5.91216 
     24.2695  
     15.3415  




```julia
S.' * randn(n)
```




    200-element Array{Float64,1}:
     23.6847 
      4.66243
      4.61544
      2.69455
     22.525  
      9.71594
     21.4546 
     14.8832 
     21.4642 
      8.61231
     29.6248 
     14.1352 
      9.6956 
      ⋮      
     11.3634 
     14.6797 
     15.4163 
      6.77538
     27.3312 
     17.3787 
     26.6867 
     10.2509 
     15.9798 
     25.6015 
     16.2479 
      8.16274



* (snparray)-matrix or matrix-(snparray) multiplications


```julia
S * randn(p, 3)
```




    100×3 Array{Float64,2}:
      -1.60062   -10.8268      5.49057 
      -9.34748     6.34819     2.06767 
       4.05002     0.325131    4.77423 
       2.9358      8.0825     -0.397737
      -5.35414   -21.3731     -1.43909 
      26.4514     -3.9009      4.35974 
       8.90361    -6.01847    11.1584  
       2.85343     0.083285    9.52995 
       2.49824    -9.35725    -3.95109 
      -0.562102  -10.0496     16.7451  
       8.27958     3.54241     2.29453 
      20.7223     -7.09494   -10.1341  
     -12.3776     -1.14019     9.97014 
       ⋮                               
      11.6659      5.81297   -35.0418  
      -6.66734   -11.7336     27.3868  
       2.32574    19.9981     22.9892  
      33.9357     -5.73126    14.894   
       0.826793   -2.0143     11.4629  
       3.18468    -6.91367     2.78527 
       5.86711     3.62248    -7.05016 
       8.42888   -16.9377     -1.87953 
      -7.8737    -10.6408      2.33838 
      -8.03682    -0.736084   10.3781  
       9.1461      8.01083    10.4097  
       8.92717    -5.15705     9.71163 




```julia
S.' * randn(n, 3)
```




    200×3 Array{Float64,2}:
     -10.6349    -19.693    -13.2774  
      -0.526305  -21.8054   -17.4738  
      -9.04862    -8.15163  -25.224   
      -3.3977    -13.7849   -16.0947  
     -10.4237    -14.0344   -14.4495  
       6.00368   -24.9057   -20.0696  
      -6.00073   -18.1709   -36.9535  
      -6.71766   -20.1876   -25.332   
     -13.1678    -22.5888   -14.3297  
       2.01221    -2.28132  -35.2057  
      -5.33609    -4.39915  -24.1266  
      -5.07446   -19.7707    -2.23233 
     -13.4401    -22.1923    -0.959362
       ⋮                              
     -12.2891     -9.20109   -4.10673 
      10.428     -13.5418   -28.1259  
     -15.5509    -23.4182   -17.6246  
      -2.16756   -33.9645    -7.70753 
     -27.8904     -4.51243  -10.2176  
     -12.3735    -29.4709   -23.3459  
     -11.7062     -9.69129  -32.2018  
      -9.16911   -29.8954   -10.581   
      -6.58973   -38.5497   -11.1196  
     -16.7095    -24.4645   -21.6888  
      -4.09502   -22.2393   -10.9904  
      -7.70413   -19.7018   -10.8812  




```julia
randn(3, n) * S
```




    3×200 Array{Any,2}:
     25.8644    11.7117   3.94609  -9.07055  …  7.97652   9.36757   9.80563
      6.44405    9.15611  5.61898  17.9318      8.86161  18.4627   23.2194 
     -7.09806  -13.1394   5.01518  -1.44749     1.2617   13.6327   -2.06252



The in-place version of matrix-(vector/matrix) multiplication functions are also available: `A_mul_B!`, `At_mul_Bt!`, `A_mul_Bt!`, `At_mul_B!`, `Ac_mul_B!`, `A_mul_Bc!`, `Ac_mul_Bc!`.


```julia
o, v = zeros(n), randn(p)
# the last optional argument is working parameter, which can preallocated to increase performance
A_mul_B!(o, S, v, similar(o))
```




    100-element Array{Float64,1}:
     26.6472 
     29.4652 
     26.6715 
     20.8541 
     24.1833 
     24.5305 
     26.5941 
     27.9289 
     35.232  
     33.375  
     21.6101 
      8.81546
     11.4446 
      ⋮      
     20.4184 
     20.4023 
     -3.50685
     25.4039 
     24.7287 
     24.0314 
     23.5559 
     25.1384 
     29.4094 
     15.9748 
     27.7337 
     25.8975 


