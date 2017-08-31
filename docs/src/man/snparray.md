
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

* `SnpArray` can be initialized from [Plink binary files](http://zzz.bwh.harvard.edu/plink/binary.shtml), say the sample data set hapmap3:


```julia
;ls -al "hapmap3.*"
```

    -rw-r--r--  1 huazhou  staff  1128171 Jun 19 14:43 hapmap3.bed
    -rw-r--r--  1 huazhou  staff   388672 Jun 19 14:43 hapmap3.bim
    -rw-r--r--  1 huazhou  staff     7136 Jun 19 14:43 hapmap3.fam
    -rw-r--r--  1 huazhou  staff   332960 Jun 19 14:43 hapmap3.map



```julia
using SnpArrays
hapmap = SnpArray("hapmap3")
```




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
hapmap = SnpArray("hapmap3"; people = 324, snps = 13928)
```




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
     (false, true)  (false, false)  (true, true)  
     (false, true)  (false, true)   (true, true)  
     (true, true)   (false, false)  (false, false)
     (true, true)   (true, true)    (false, true) 
     (true, true)   (false, true)   (false, true) 



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




    (true, true)



`randgeno(n, maf, minor_allele)` generates a vector of random genotypes according to a common minor allele frequency `maf` and the minor allele.


```julia
randgeno(10, 0.25, true)
```




    10-element SnpArrays.SnpArray{1}:
     (false, true) 
     (true, true)  
     (true, true)  
     (false, false)
     (true, true)  
     (true, true)  
     (true, true)  
     (false, true) 
     (false, false)
     (false, true) 



`randgeno(m, n, maf, minor_allele)` generates a random $m$-by-$n$ `SnpArray` according to a vector of minor allele frequencies `maf` and a minor allele indicator vector. The lengths of both vectors should be `n`.


```julia
# this is a random replicate of the hapmap data
randgeno(size(hapmap), maf, minor_allele)
```




    324×13928 SnpArrays.SnpArray{2}:
     (true, true)  (false, true)  (true, true)    …  (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (false, true)  (true, true)
     (true, true)  (false, true)  (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)    …  (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)       (false, true)  (true, true)
     (true, true)  (false, true)  (false, true)      (true, true)   (true, true)
     (true, true)  (false, true)  (true, true)       (true, true)   (true, true)
     (true, true)  (false, true)  (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)    …  (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)
     ⋮                                            ⋱                             
     (true, true)  (false, true)  (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)       (false, true)  (true, true)
     (true, true)  (true, true)   (false, true)   …  (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)     (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)   …  (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)



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
     0.571346   0.0465045  0.0204832  …  0.0650646  0.0716494  0.0650761
     0.0465045  0.544773   0.028394      0.0513636  0.0449712  0.0631397
     0.0204832  0.028394   0.519855      0.0457716  0.0301285  0.0364614
     0.0481785  0.0370097  0.0285135     0.0599191  0.0657765  0.0599049
     0.0527814  0.0431142  0.0255107     0.0722031  0.0583019  0.0659732
     0.0446293  0.031857   0.039093   …  0.0703164  0.0569839  0.065417 
     0.0398707  0.0225306  0.0122048     0.044638   0.0383134  0.0375679
     0.0415416  0.0383736  0.021525      0.0576784  0.0553812  0.0662787
     0.030548   0.0313828  0.0164342     0.0336193  0.0465423  0.0380742
     0.0392845  0.0417486  0.0257594     0.0666458  0.0577377  0.0490077
     0.0479005  0.0456751  0.0231647  …  0.0584296  0.06287    0.0610135
     0.061249   0.0390106  0.0371638     0.0660586  0.0582874  0.0701743
     0.0360434  0.0440139  0.0266148     0.0574341  0.070172   0.063768 
     ⋮                                ⋱                                 
     0.0582425  0.0552193  0.0394548     0.0863796  0.0865721  0.0968307
     0.067412   0.0522998  0.033598      0.0918259  0.0802092  0.0854516
     0.0680476  0.0560847  0.0388254     0.0899582  0.08414    0.0854648
     0.0656848  0.0569648  0.039208   …  0.0877947  0.0764744  0.0892505
     0.0626458  0.0555299  0.0368724     0.0785841  0.0735793  0.0814078
     0.0649787  0.0605798  0.0406489     0.0889447  0.0902569  0.0812508
     0.0624551  0.0603942  0.0388867     0.0825189  0.0849712  0.0817198
     0.0621368  0.0631124  0.0326863     0.092926   0.0901677  0.0864766
     0.0660474  0.0631126  0.032413   …  0.0975443  0.0798605  0.0866241
     0.0650646  0.0513636  0.0457716     0.603213   0.0820706  0.0905677
     0.0716494  0.0449712  0.0301285     0.0820706  0.595356   0.0824803
     0.0650761  0.0631397  0.0364614     0.0905677  0.0824803  0.583244 




```julia
# GRM using all SNPs with MAF ≥ 0.05
grm(hapmap; maf_threshold = 0.05)
```




    324×324 Array{Float64,2}:
     0.571452   0.0462904  0.0199165  …  0.0651794  0.0717623  0.0651457
     0.0462904  0.544914   0.0283713     0.0511696  0.0456596  0.0630129
     0.0199165  0.0283713  0.520022      0.0460785  0.030228   0.0359554
     0.0485764  0.0370786  0.0292213     0.05997    0.0658113  0.060026 
     0.0529626  0.0429542  0.0254827     0.072063   0.0589913  0.0658106
     0.04432    0.0314165  0.0385681  …  0.0707138  0.0562916  0.0652219
     0.0397146  0.0225996  0.0116645     0.0443248  0.0382372  0.0371071
     0.0410736  0.0386706  0.0215219     0.0577337  0.0546118  0.0661473
     0.0301126  0.0308084  0.0167159     0.0335706  0.0458171  0.0378848
     0.0395036  0.0420549  0.0263567     0.0665525  0.0577466  0.0489306
     0.048092   0.0458081  0.022875   …  0.0585571  0.0639564  0.0613725
     0.0607261  0.0389591  0.0369648     0.066402   0.0582424  0.0700461
     0.036022   0.0441741  0.0264965     0.0574379  0.0701088  0.0632619
     ⋮                                ⋱                                 
     0.0577144  0.0552319  0.0385214     0.0862471  0.086786   0.0971796
     0.0672011  0.0521327  0.0330645     0.0917182  0.0804224  0.0859455
     0.0682613  0.0557235  0.0386046     0.0897848  0.0842798  0.0852609
     0.0652514  0.057126   0.0391473  …  0.0874249  0.0764886  0.089188 
     0.0626828  0.0563047  0.0361525     0.0787534  0.0731322  0.0818519
     0.0652242  0.0609512  0.0402041     0.0891027  0.0902155  0.0810412
     0.0624822  0.0607759  0.0382927     0.0831698  0.0854469  0.0817415
     0.0622411  0.0631895  0.0334202     0.0930534  0.0904701  0.0866312
     0.0657874  0.063138   0.0328947  …  0.0972887  0.0794073  0.086289 
     0.0651794  0.0511696  0.0460785     0.603223   0.0820763  0.0904615
     0.0717623  0.0456596  0.030228      0.0820763  0.595456   0.0823627
     0.0651457  0.0630129  0.0359554     0.0904615  0.0823627  0.583026 




```julia
# GRM using every other SNP, with maf ≥ 0.01
grm(view(hapmap, :, 1:2:snps))
```




    324×324 Array{Float64,2}:
     0.559196   0.0430886  0.0275476  …  0.0675327  0.0740519  0.068678 
     0.0430886  0.558878   0.0286428     0.05786    0.045359   0.0561938
     0.0275476  0.0286428  0.513396      0.038853   0.0392725  0.0475248
     0.0457033  0.0459737  0.0272105     0.0514827  0.0619385  0.0564087
     0.0521454  0.0481397  0.0257591     0.0672912  0.0569231  0.0659968
     0.0522916  0.041013   0.0398454  …  0.0773389  0.0620592  0.0539301
     0.0397736  0.0267719  0.016085      0.0471858  0.0346203  0.03449  
     0.047408   0.0376607  0.0259769     0.0589596  0.0556453  0.0687727
     0.0269957  0.0228234  0.0201158     0.0332562  0.0440439  0.0361991
     0.0321618  0.0401614  0.0232686     0.0601785  0.0409768  0.0504197
     0.0494131  0.050353   0.0212866  …  0.0629188  0.0651026  0.053264 
     0.062188   0.0503319  0.0407011     0.0663986  0.0647543  0.0607589
     0.0253532  0.0426208  0.0228778     0.0549272  0.0689409  0.0634824
     ⋮                                ⋱                                 
     0.061482   0.0516746  0.0410709     0.0870478  0.0848473  0.102271 
     0.0664257  0.0600823  0.03027       0.0891805  0.0758972  0.0808276
     0.0692706  0.0538586  0.0393204     0.0847024  0.081855   0.096149 
     0.0651635  0.0578905  0.0415709  …  0.0936544  0.0772461  0.0823608
     0.063133   0.0625333  0.042042      0.0730644  0.0690764  0.0872862
     0.0663681  0.0638972  0.0420809     0.077098   0.0829472  0.0854734
     0.0650961  0.0637289  0.0307002     0.0862602  0.092594   0.0847355
     0.059969   0.0619758  0.0343407     0.0949796  0.0913884  0.0837603
     0.076881   0.0664363  0.0395719  …  0.101928   0.0746607  0.0900927
     0.0675327  0.05786    0.038853      0.573711   0.0699417  0.095876 
     0.0740519  0.045359   0.0392725     0.0699417  0.57733    0.074765 
     0.068678   0.0561938  0.0475248     0.095876   0.074765   0.588363 




```julia
# MoM using all SNPs with MAF ≥ 0.01
grm(hapmap; method = :MoM)
```




    324×324 Array{Float64,2}:
     0.539203   0.0350323  0.0024137   …  0.0539417  0.0638691  0.0509871
     0.0350323  0.517812   0.0136411      0.0422415  0.0397597  0.0507507
     0.0024137  0.0136411  0.499966       0.0327868  0.0220322  0.0196685
     0.0436597  0.0291232  0.0239231      0.0526417  0.0695419  0.0495689
     0.0449597  0.0331414  0.0152957      0.0654055  0.0568963  0.0567781
     0.0323141  0.0202594  0.0260504   …  0.0600872  0.0498053  0.0493325
     0.0251049  0.0113956  0.00371372     0.0311323  0.0308959  0.0247504
     0.0255777  0.0280595  0.00891379     0.0385778  0.0462598  0.0463779
     0.0216776  0.0256958  0.0136411      0.0301868  0.048978   0.0304232
     0.0210867  0.0247504  0.0104502      0.0524053  0.0509871  0.0319596
     0.0359778  0.0331414  0.00595921  …  0.0499234  0.0563053  0.0505144
     0.0485052  0.0324323  0.0269959      0.0559508  0.0555962  0.0521689
     0.0314868  0.0421233  0.0247504      0.0589054  0.0708419  0.0579599
     ⋮                                 ⋱                                 
     0.0421233  0.0417688  0.025814       0.0713146  0.0793511  0.0821875
     0.0624509  0.0481507  0.0269959      0.0916422  0.0793511  0.083133 
     0.0550053  0.0455507  0.0255777      0.0806511  0.0817148  0.0726147
     0.0453143  0.0428324  0.0240413   …  0.072851   0.0625691  0.0704874
     0.056069   0.0525235  0.0292413      0.0771056  0.0704874  0.074151 
     0.0529962  0.0550053  0.0314868      0.0788784  0.0869149  0.0726147
     0.0526417  0.0508689  0.0297141      0.0780511  0.0834875  0.0732056
     0.0487416  0.0524053  0.0212049      0.0821875  0.0843148  0.0726147
     0.0503962  0.0477961  0.0213231   …  0.0880967  0.0737965  0.0691873
     0.0539417  0.0422415  0.0327868      0.561421   0.081242   0.076042 
     0.0638691  0.0397597  0.0220322      0.081242   0.567094   0.0693055
     0.0509871  0.0507507  0.0196685      0.076042   0.0693055  0.533766 



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


