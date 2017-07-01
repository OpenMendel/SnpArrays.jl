
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

    -rw-r--r--  1 huazhou  staff  1128171 Jun 30 07:48 hapmap3.bed
    -rw-r--r--  1 huazhou  staff   388672 Jun 30 07:48 hapmap3.bim
    -rw-r--r--  1 huazhou  staff     7136 Jun 30 07:48 hapmap3.fam
    -rw-r--r--  1 huazhou  staff   332960 Jun 30 07:48 hapmap3.map



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
     (false, false)  (false, false)  (false, true) 
     (true, true)    (true, true)    (false, true) 
     (false, true)   (false, false)  (false, false)
     (false, false)  (false, true)   (false, true) 
     (true, true)    (false, false)  (false, true) 



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




    (false, true)



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
     (false, true)
     (false, true)
     (false, true)
     (false, true)
     (true, true) 
     (false, true)
     (true, true) 



`randgeno(m, n, maf, minor_allele)` generates a random $m$-by-$n$ `SnpArray` according to a vector of minor allele frequencies `maf` and a minor allele indicator vector. The lengths of both vectors should be `n`.


```julia
# this is a random replicate of the hapmap data
randgeno(size(hapmap), maf, minor_allele)
```




    324×13928 SnpArrays.SnpArray{2}:
     (true, true)  (true, true)   (true, true)    …  (true, true)   (true, true)
     (true, true)  (false, true)  (true, true)       (false, true)  (true, true)
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)
     (true, true)  (false, true)  (true, true)       (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)    …  (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (false, true)  (true, true)       (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)   …  (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (false, true)  (true, true)       (true, true)   (true, true)
     ⋮                                            ⋱                             
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (false, true)  (true, true)       (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)    …  (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (false, true)  (true, true)
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)
     (true, true)  (true, true)   (false, false)     (true, true)   (true, true)
     (true, true)  (true, true)   (false, true)      (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)    …  (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)
     (true, true)  (true, true)   (true, true)       (true, true)   (true, true)
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
    



```julia
# GRM using all SNPs
grm(hapmap)
```




    324×324 Array{Float64,2}:
     0.566466   0.0444359  0.0190432  …  0.0626215  0.0687383  0.0622614
     0.0444359  0.530419   0.0309971     0.0495598  0.0432578  0.0603729
     0.0190432  0.0309971  0.511297      0.0446246  0.0294237  0.035185 
     0.0462825  0.0352154  0.0275386     0.0575497  0.0632928  0.0572074
     0.0508401  0.041172   0.0238817     0.0693205  0.0560045  0.06329  
     0.042911   0.0304459  0.037073   …  0.0680288  0.0542302  0.0625457
     0.0379732  0.0212719  0.012149      0.0427044  0.0365571  0.0355025
     0.0396982  0.0372525  0.0210437     0.0557044  0.05327    0.0632351
     0.0288591  0.0298873  0.0163501     0.032368   0.044781   0.0364853
     0.0373909  0.0410606  0.0251445     0.0641125  0.0553558  0.0466605
     0.0458546  0.0439146  0.0224671  …  0.0562275  0.0643399  0.0586755
     0.0580315  0.0368683  0.0352559     0.0636779  0.0561118  0.0689767
     0.0347037  0.0423452  0.0259498     0.0550224  0.0675027  0.0608943
     ⋮                                ⋱                                 
     0.0633259  0.0531724  0.0371187     0.0829597  0.0831308  0.0931435
     0.0645506  0.0498678  0.0329406     0.0883068  0.0772216  0.0853762
     0.0653716  0.0534623  0.0372739     0.0863184  0.0812877  0.0821439
     0.0627039  0.0548826  0.0383003  …  0.0846339  0.0734946  0.0857409
     0.0602445  0.0536981  0.0355785     0.0755747  0.0704755  0.0783867
     0.062317   0.0589174  0.0383417     0.0850735  0.0887094  0.0795435
     0.0602235  0.0586656  0.0367576     0.0793016  0.0808723  0.0779107
     0.0594668  0.0606838  0.0316561     0.0895881  0.0864562  0.0831714
     0.0635207  0.0607507  0.031338   …  0.0937335  0.076364   0.0831758
     0.0626215  0.0495598  0.0446246     0.606903   0.0789186  0.0868736
     0.0687383  0.0432578  0.0294237     0.0789186  0.583857   0.0792285
     0.0622614  0.0603729  0.035185      0.0868736  0.0792285  0.5758   




```julia
# GRM using every other SNP
grm(view(hapmap, :, 1:2:snps))
```




    324×324 Array{Float64,2}:
     0.555683   0.0415085  0.0264839  …  0.0651673  0.0714261  0.0656155
     0.0415085  0.545431   0.0354667     0.0562859  0.0438491  0.0536016
     0.0264839  0.0354667  0.500566      0.0375265  0.0371542  0.0453827
     0.0434207  0.0446201  0.0253641     0.0492325  0.059049   0.0542097
     0.0499442  0.046307   0.0248879     0.0655328  0.0548686  0.062911 
     0.0502347  0.0392831  0.038141   …  0.0739908  0.0596454  0.0509341
     0.03778    0.02625    0.0156034     0.0449416  0.0331082  0.0330955
     0.045428   0.0375521  0.0252809     0.0569317  0.054053   0.0663885
     0.0253056  0.023008   0.0177888     0.0326947  0.0417944  0.0346796
     0.0308492  0.0385291  0.0232155     0.057392   0.0396136  0.0486733
     0.0473886  0.0489292  0.0200301  …  0.061013   0.0623889  0.0514638
     0.0606601  0.047154   0.0381714     0.0635991  0.0628481  0.062388 
     0.0239943  0.041524   0.0214056     0.0527117  0.0659221  0.0610669
     ⋮                                ⋱                                 
     0.0585219  0.0484891  0.0385751     0.0841442  0.0819137  0.0982057
     0.0640031  0.0582749  0.0286075     0.0857168  0.073329   0.0839684
     0.066651   0.0521818  0.0373679     0.0813776  0.0801579  0.0922047
     0.0626328  0.0566101  0.0382988  …  0.0902857  0.073534   0.0799175
     0.0613747  0.0595828  0.0398474     0.0700391  0.0686345  0.0845796
     0.063787   0.0615761  0.0397324     0.0731029  0.0793144  0.085481 
     0.0617635  0.0614787  0.0298845     0.0825804  0.0893492  0.0811897
     0.0581062  0.0591602  0.0323554     0.0917713  0.0875897  0.0803352
     0.0737464  0.0643129  0.0385646  …  0.098534   0.0733618  0.0866629
     0.0651673  0.0562859  0.0375265     0.600007   0.0672816  0.0920241
     0.0714261  0.0438491  0.0371542     0.0672816  0.562885   0.0721028
     0.0656155  0.0536016  0.0453827     0.0920241  0.0721028  0.576308 




```julia
# MoM using all SNPs
grm(hapmap; method = :MoM)
```




    324×324 Array{Float64,2}:
     0.53945     0.0347339  0.00344015  …  0.0535102  0.0631936  0.0506761
     0.0347339   0.517957   0.0150129      0.0420555  0.0395756  0.0497313
     0.00344015  0.0150129  0.49989        0.033435   0.0223345  0.0206813
     0.0428821   0.0289475  0.0239878      0.0519751  0.0690981  0.0490228
     0.0448897   0.0333169  0.0158396      0.0649649  0.0564625  0.0555177
     0.0320179   0.0215079  0.0267038   …  0.0598871  0.0502037  0.0494952
     0.0248144   0.0112341  0.00379442     0.0313093  0.0306008  0.0245782
     0.0262315   0.0290656  0.0105255      0.0389852  0.0464248  0.0467791
     0.0209174   0.0261134  0.0140682      0.0294199  0.0486685  0.0304827
     0.0212717   0.0259953  0.0117064      0.053156   0.0512665  0.0322541
     0.0356787   0.0336711  0.00651048  …  0.0494952  0.0559901  0.0500856
     0.04796     0.0326083  0.0277666      0.0552816  0.0552816  0.0523293
     0.0309551   0.0421736  0.0255229      0.0582338  0.070279   0.0575253
     ⋮                                  ⋱                                 
     0.0419374   0.0419374  0.0261134      0.0712237  0.0792538  0.0814975
     0.0621308   0.0481962  0.0281209      0.091417   0.0790176  0.0832689
     0.054573    0.0455982  0.0267038      0.0804347  0.0810251  0.0727589
     0.0457163   0.0427641  0.0248144   …  0.0726408  0.0627212  0.0699247
     0.0559901   0.0533921  0.0301284      0.0777186  0.070279   0.0744121
     0.0516208   0.0544549  0.0313093      0.0790176  0.0857487  0.0713418
     0.0513846   0.0506761  0.0300104      0.0773644  0.0830327  0.0738217
     0.0490228   0.0526836  0.0219802      0.0816156  0.0839774  0.0718141
     0.0504399   0.0481962  0.0226888   …  0.0883467  0.0738217  0.0692162
     0.0535102   0.0420555  0.033435       0.561533   0.080789   0.0764196
     0.0631936   0.0395756  0.0223345      0.080789   0.56661    0.06898  
     0.0506761   0.0497313  0.0206813      0.0764196  0.06898    0.533545 



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




    ([-38.7231 1.2983 -7.00541; -32.6096 1.21052 -3.3232; … ; -48.9263 2.06102 2.17374; -48.8627 -0.274894 6.49518], [-1.64162e-18 7.41502e-19 -5.96439e-18; 0.00143962 0.0042375 -0.00311816; … ; 0.00313326 0.00427486 -0.0152038; -9.09523e-5 0.00287777 0.0037855], [1841.4, 225.324, 70.7084])



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
      225.323 
       70.7085



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

    [1m[33mWARNING: [39m[22m[33mArray{T}(::Type{T}, m::Int) is deprecated, use Array{T}(m) instead.[39m
    Stacktrace:
     [1] [1mdepwarn[22m[22m[1m([22m[22m::String, ::Symbol[1m)[22m[22m at [1m./deprecated.jl:70[22m[22m
     [2] [1mArray[22m[22m[1m([22m[22m::Type{Float32}, ::Int64[1m)[22m[22m at [1m./deprecated.jl:57[22m[22m
     [3] [1mA_mul_B![22m[22m[1m([22m[22m::SubArray{Float32,1,Array{Float32,1},Tuple{UnitRange{Int64}},true}, ::LinearMaps.CompositeMap{Float32,Tuple{LinearMaps.FunctionMap{Float32,SnpArrays.##10#12{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}},SnpArrays.##11#13{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}}},LinearMaps.TransposeMap{Float32,LinearMaps.FunctionMap{Float32,SnpArrays.##10#12{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}},SnpArrays.##11#13{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}}}}}}, ::SubArray{Float32,1,Array{Float32,1},Tuple{UnitRange{Int64}},true}[1m)[22m[22m at [1m/Users/huazhou/.julia/v0.6/LinearMaps/src/composition.jl:88[22m[22m
     [4] [1maupd_wrapper[22m[22m[1m([22m[22m::Type{T} where T, ::Base.LinAlg.#matvecA!#114{LinearMaps.CompositeMap{Float32,Tuple{LinearMaps.FunctionMap{Float32,SnpArrays.##10#12{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}},SnpArrays.##11#13{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}}},LinearMaps.TransposeMap{Float32,LinearMaps.FunctionMap{Float32,SnpArrays.##10#12{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}},SnpArrays.##11#13{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}}}}}}}, ::Base.LinAlg.##108#115, ::Base.LinAlg.##109#116, ::Int64, ::Bool, ::Bool, ::String, ::Int64, ::Int64, ::String, ::Float64, ::Int64, ::Int64, ::Array{Float32,1}[1m)[22m[22m at [1m./linalg/arpack.jl:59[22m[22m
     [5] [1m#_eigs#107[22m[22m[1m([22m[22m::Int64, ::Int64, ::Symbol, ::Float64, ::Int64, ::Void, ::Array{Float32,1}, ::Bool, ::Base.LinAlg.#_eigs, ::LinearMaps.CompositeMap{Float32,Tuple{LinearMaps.FunctionMap{Float32,SnpArrays.##10#12{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}},SnpArrays.##11#13{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}}},LinearMaps.TransposeMap{Float32,LinearMaps.FunctionMap{Float32,SnpArrays.##10#12{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}},SnpArrays.##11#13{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}}}}}}, ::UniformScaling{Int64}[1m)[22m[22m at [1m./linalg/arnoldi.jl:285[22m[22m
     [6] [1m(::Base.LinAlg.#kw##_eigs)[22m[22m[1m([22m[22m::Array{Any,1}, ::Base.LinAlg.#_eigs, ::LinearMaps.CompositeMap{Float32,Tuple{LinearMaps.FunctionMap{Float32,SnpArrays.##10#12{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}},SnpArrays.##11#13{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}}},LinearMaps.TransposeMap{Float32,LinearMaps.FunctionMap{Float32,SnpArrays.##10#12{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}},SnpArrays.##11#13{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}}}}}}, ::UniformScaling{Int64}[1m)[22m[22m at [1m./<missing>:0[22m[22m
     [7] [1m#eigs#106[22m[22m[1m([22m[22m::Array{Any,1}, ::Function, ::LinearMaps.CompositeMap{Float32,Tuple{LinearMaps.FunctionMap{Float32,SnpArrays.##10#12{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}},SnpArrays.##11#13{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}}},LinearMaps.TransposeMap{Float32,LinearMaps.FunctionMap{Float32,SnpArrays.##10#12{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}},SnpArrays.##11#13{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}}}}}}, ::UniformScaling{Int64}[1m)[22m[22m at [1m./linalg/arnoldi.jl:170[22m[22m
     [8] [1m(::Base.LinAlg.#kw##eigs)[22m[22m[1m([22m[22m::Array{Any,1}, ::Base.LinAlg.#eigs, ::LinearMaps.CompositeMap{Float32,Tuple{LinearMaps.FunctionMap{Float32,SnpArrays.##10#12{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}},SnpArrays.##11#13{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}}},LinearMaps.TransposeMap{Float32,LinearMaps.FunctionMap{Float32,SnpArrays.##10#12{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}},SnpArrays.##11#13{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}}}}}}, ::UniformScaling{Int64}[1m)[22m[22m at [1m./<missing>:0[22m[22m
     [9] [1m#eigs#99[22m[22m at [1m./linalg/arnoldi.jl:90[22m[22m [inlined]
     [10] [1m(::Base.LinAlg.#kw##eigs)[22m[22m[1m([22m[22m::Array{Any,1}, ::Base.LinAlg.#eigs, ::LinearMaps.CompositeMap{Float32,Tuple{LinearMaps.FunctionMap{Float32,SnpArrays.##10#12{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}},SnpArrays.##11#13{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}}},LinearMaps.TransposeMap{Float32,LinearMaps.FunctionMap{Float32,SnpArrays.##10#12{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}},SnpArrays.##11#13{SparseMatrixCSC{Float32,UInt32},Array{Float32,1},Array{Float32,1}}}}}}[1m)[22m[22m at [1m./<missing>:0[22m[22m
     [11] [1mpca_sp[22m[22m[1m([22m[22m::SnpArrays.SnpArray{2}, ::Int64, ::Type{SparseMatrixCSC{Float32,UInt32}}[1m)[22m[22m at [1m/Users/huazhou/.julia/v0.6/SnpArrays/src/SnpArrays.jl:652[22m[22m
     [12] [1minclude_string[22m[22m[1m([22m[22m::String, ::String[1m)[22m[22m at [1m./loading.jl:515[22m[22m
     [13] [1mexecute_request[22m[22m[1m([22m[22m::ZMQ.Socket, ::IJulia.Msg[1m)[22m[22m at [1m/Users/huazhou/.julia/v0.6/IJulia/src/execute_request.jl:160[22m[22m
     [14] [1meventloop[22m[22m[1m([22m[22m::ZMQ.Socket[1m)[22m[22m at [1m/Users/huazhou/.julia/v0.6/IJulia/src/eventloop.jl:8[22m[22m
     [15] [1m(::IJulia.##11#14)[22m[22m[1m([22m[22m[1m)[22m[22m at [1m./task.jl:335[22m[22m
    while loading In[94], in expression starting on line 3





    3-element Array{Float32,1}:
     1841.39  
      225.324 
       70.7083


