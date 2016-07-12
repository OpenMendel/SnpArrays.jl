
# SnpArray

`SnpArray` is an array of `Tuple{Bool,Bool}` and adopts the same coding as the [Plink binary format](http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml). If `A1` and `A2` are the two alleles, the coding rule is  

| Genotype | SnpArray |  
|---|---|---|  
| A1,A1 | (false,false) |  
| A1,A2 | (false,true) |  
| A2,A2 | (true,true) |  
| missing | (true,false) |  

The code `(true,false)` is reserved for missing genotype. Otherwise, the bit `true` represents one copy of allele `A2`. In a two-dimensional `SnpArray`, each row is a person and each column is a SNP.

For complete genotype data, for example, after imputation, consider using the [`HaplotypeArray`](https://github.com/OpenMendel/SnpArrays.jl/blob/master/docs/haplotypearray.ipynb) type.

## Constructor

There are various ways to initialize a `SnpArray`.  

* `SnpArray` can be initialized from [Plink binary files](http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml), say the sample data set hapmap3:


```julia
;ls -al hapmap3.*
```

    -rw-r--r--  1 hzhou3  staff  1128171 Jul 11 15:43 hapmap3.bed
    -rw-r--r--  1 hzhou3  staff   388672 Jul  7 18:09 hapmap3.bim
    -rw-r--r--  1 hzhou3  staff     7136 Jul  7 18:09 hapmap3.fam
    -rw-r--r--  1 hzhou3  staff   332960 Jul  7 18:09 hapmap3.map



```julia
using SnpArrays
hapmap = SnpArray("hapmap3")
```




    324x13928 SnpArrays.SnpArray{2}:
     (true,true)  (true,true)   (false,false)  …  (true,true)   (true,true)
     (true,true)  (false,true)  (false,true)      (false,true)  (true,true)
     (true,true)  (true,true)   (false,true)      (true,true)   (true,true)
     (true,true)  (true,true)   (false,true)      (true,true)   (true,true)
     (true,true)  (true,true)   (false,true)      (true,true)   (true,true)
     (true,true)  (false,true)  (true,true)    …  (true,true)   (true,true)
     (true,true)  (true,true)   (true,true)       (true,true)   (true,true)
     (true,true)  (true,true)   (false,false)     (true,true)   (true,true)
     (true,true)  (true,true)   (false,true)      (true,true)   (true,true)
     (true,true)  (true,true)   (false,true)      (true,true)   (true,true)
     (true,true)  (true,true)   (false,true)   …  (true,true)   (true,true)
     (true,true)  (true,true)   (true,true)       (true,true)   (true,true)
     (true,true)  (true,true)   (false,false)     (true,true)   (true,true)
     ⋮                                         ⋱                           
     (true,true)  (true,true)   (false,true)      (true,true)   (true,true)
     (true,true)  (true,true)   (false,false)     (true,true)   (true,true)
     (true,true)  (true,true)   (false,false)     (true,true)   (true,true)
     (true,true)  (true,true)   (false,false)  …  (true,true)   (true,true)
     (true,true)  (true,true)   (false,true)      (true,true)   (true,true)
     (true,true)  (true,true)   (true,true)       (true,true)   (true,true)
     (true,true)  (true,true)   (false,true)      (true,true)   (true,true)
     (true,true)  (true,true)   (false,true)      (true,true)   (true,true)
     (true,true)  (true,true)   (false,false)  …  (true,true)   (true,true)
     (true,true)  (true,true)   (false,true)      (true,true)   (true,true)
     (true,true)  (true,true)   (false,false)     (true,true)   (true,true)
     (true,true)  (true,true)   (false,false)     (true,true)   (true,true)



By default, the constructor figures out the number of individuals and SNPs from the `.bim` and `.fam` files.


```julia
# rows are people; columns are SNPs
people, snps = size(hapmap)
```




    (324,13928)



Alternatively, users can supply keyword arguments `people` and `snps` to the constructor. In this case only the `.bed` file needs to be present.


```julia
hapmap = SnpArray("hapmap3"; people=324, snps=13928)
```




    324x13928 SnpArrays.SnpArray{2}:
     (true,true)  (true,true)   (false,false)  …  (true,true)   (true,true)
     (true,true)  (false,true)  (false,true)      (false,true)  (true,true)
     (true,true)  (true,true)   (false,true)      (true,true)   (true,true)
     (true,true)  (true,true)   (false,true)      (true,true)   (true,true)
     (true,true)  (true,true)   (false,true)      (true,true)   (true,true)
     (true,true)  (false,true)  (true,true)    …  (true,true)   (true,true)
     (true,true)  (true,true)   (true,true)       (true,true)   (true,true)
     (true,true)  (true,true)   (false,false)     (true,true)   (true,true)
     (true,true)  (true,true)   (false,true)      (true,true)   (true,true)
     (true,true)  (true,true)   (false,true)      (true,true)   (true,true)
     (true,true)  (true,true)   (false,true)   …  (true,true)   (true,true)
     (true,true)  (true,true)   (true,true)       (true,true)   (true,true)
     (true,true)  (true,true)   (false,false)     (true,true)   (true,true)
     ⋮                                         ⋱                           
     (true,true)  (true,true)   (false,true)      (true,true)   (true,true)
     (true,true)  (true,true)   (false,false)     (true,true)   (true,true)
     (true,true)  (true,true)   (false,false)     (true,true)   (true,true)
     (true,true)  (true,true)   (false,false)  …  (true,true)   (true,true)
     (true,true)  (true,true)   (false,true)      (true,true)   (true,true)
     (true,true)  (true,true)   (true,true)       (true,true)   (true,true)
     (true,true)  (true,true)   (false,true)      (true,true)   (true,true)
     (true,true)  (true,true)   (false,true)      (true,true)   (true,true)
     (true,true)  (true,true)   (false,false)  …  (true,true)   (true,true)
     (true,true)  (true,true)   (false,true)      (true,true)   (true,true)
     (true,true)  (true,true)   (false,false)     (true,true)   (true,true)
     (true,true)  (true,true)   (false,false)     (true,true)   (true,true)



Internally `SnpArray` stores data as `BitArray`s and consumes approximately the same amount of memory as the Plink `bed` file size.


```julia
# memory usage, bed file size
Base.summarysize(hapmap), filesize("hapmap3.bed")
```




    (1128256,1128171)



* `SnpArray` can be initialized from a matrix of A1 allele counts.


```julia
SnpArray(rand(0:2, 5, 3))
```




    5x3 SnpArrays.SnpArray{2}:
     (false,false)  (false,true)  (false,false)
     (true,true)    (false,true)  (true,true)  
     (false,true)   (false,true)  (false,false)
     (true,true)    (false,true)  (false,false)
     (false,false)  (true,true)   (false,true) 



* `SnpArray(m, n)` generates an m by n `SnpArray` of all A1 alleles.


```julia
s = SnpArray(5, 3)
```




    5x3 SnpArrays.SnpArray{2}:
     (false,false)  (false,false)  (false,false)
     (false,false)  (false,false)  (false,false)
     (false,false)  (false,false)  (false,false)
     (false,false)  (false,false)  (false,false)
     (false,false)  (false,false)  (false,false)



## Summary statistics

`summarize` function computes the following summary statistics of a `SnpArray`  
* `maf`: minor allele frequencies, taking into account of missingness.  
* `minor_allele`: a `BitVector` indicating the minor allele for each SNP. `minor_allele[j]==true` means A1 is the minor allele for SNP j; `minor_allele[j]==false` means A2 is the minor allele for SNP j.  
* `missings_by_snp`: number of missing genotypes for each snp.  
* `missings_by_person`: number of missing genotypes for each person.


```julia
maf, minor_allele, missings_by_snp, missings_by_person = summarize(hapmap)
# minor allele frequencies
maf'
```




    1x13928 Array{Float64,2}:
     0.0  0.0776398  0.324074  0.191589  …  0.00154321  0.0417957  0.00617284




```julia
# total number of missing genotypes
sum(missings_by_snp), sum(missings_by_person)
```




    (11894,11894)




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
_, _, missings_by_snp_filtered, missings_by_person_filtered = summarize(sub(hapmap, person_idx, snp_idx));
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




    (false,true)



`randgeno(maf, minor_allele)` generates a random genotype according to minor allele frequency `maf` and whether the minor allele is A1 (`minor_allele==true`) or A2 (`minor_allele==false`).


```julia
randgeno(0.25, true)
```




    (true,true)



`randgeno(n, maf, minor_allele)` generates a vector of random genotypes according to a common minor allele frequency `maf` and the minor allele.


```julia
randgeno(10, 0.25, true)
```




    10-element SnpArrays.SnpArray{1}:
     (false,true)
     (false,true)
     (true,true) 
     (true,true) 
     (false,true)
     (false,true)
     (false,true)
     (false,true)
     (false,true)
     (true,true) 



`randgeno(m, n, maf, minor_allele)` generates a random $m$-by-$n$ `SnpArray` according to a vector of minor allele frequencies `maf` and a minor allele indicator vector. The lengths of both vectors should be `n`.


```julia
# this is a random replicate of the hapmap data
randgeno(size(hapmap), maf, minor_allele)
```




    324x13928 SnpArrays.SnpArray{2}:
     (true,true)  (true,true)    (true,true)   …  (false,true)  (true,true)
     (true,true)  (false,true)   (true,true)      (true,true)   (true,true)
     (true,true)  (true,true)    (true,true)      (true,true)   (true,true)
     (true,true)  (true,true)    (true,true)      (true,true)   (true,true)
     (true,true)  (false,true)   (false,true)     (true,true)   (true,true)
     (true,true)  (true,true)    (true,true)   …  (true,true)   (true,true)
     (true,true)  (true,true)    (false,true)     (true,true)   (true,true)
     (true,true)  (false,false)  (true,true)      (true,true)   (true,true)
     (true,true)  (true,true)    (true,true)      (true,true)   (true,true)
     (true,true)  (false,true)   (false,true)     (true,true)   (true,true)
     (true,true)  (true,true)    (true,true)   …  (true,true)   (true,true)
     (true,true)  (true,true)    (false,true)     (true,true)   (true,true)
     (true,true)  (true,true)    (false,true)     (true,true)   (true,true)
     ⋮                                         ⋱                           
     (true,true)  (true,true)    (false,true)     (true,true)   (true,true)
     (true,true)  (true,true)    (true,true)      (true,true)   (true,true)
     (true,true)  (true,true)    (true,true)      (true,true)   (true,true)
     (true,true)  (true,true)    (true,true)   …  (false,true)  (true,true)
     (true,true)  (false,true)   (true,true)      (true,true)   (true,true)
     (true,true)  (true,true)    (false,true)     (true,true)   (true,true)
     (true,true)  (true,true)    (false,true)     (false,true)  (true,true)
     (true,true)  (true,true)    (true,true)      (true,true)   (true,true)
     (true,true)  (true,true)    (false,true)  …  (true,true)   (true,true)
     (true,true)  (false,true)   (true,true)      (true,true)   (true,true)
     (true,true)  (true,true)    (true,true)      (true,true)   (true,true)
     (true,true)  (true,true)    (true,true)      (true,true)   (true,true)



## Subsetting

Subsetting a `SnpArray` works the same way as subsetting any other arrays.


```julia
# genotypes of the 1st person
hapmap[1, :]
```




    1x13928 SnpArrays.SnpArray{2}:
     (true,true)  (true,true)  (false,false)  …  (true,true)  (true,true)




```julia
# genotypes of the 5th SNP
hapmap[:, 5]
```




    324-element SnpArrays.SnpArray{1}:
     (true,true)  
     (true,true)  
     (false,true) 
     (false,true) 
     (true,true)  
     (false,false)
     (false,false)
     (true,true)  
     (true,true)  
     (true,true)  
     (true,true)  
     (true,true)  
     (false,true) 
     ⋮            
     (false,false)
     (true,true)  
     (false,true) 
     (true,true)  
     (true,true)  
     (true,true)  
     (true,true)  
     (true,true)  
     (false,true) 
     (true,true)  
     (true,true)  
     (true,true)  




```julia
# subsetting both persons and SNPs
hapmap[1:5, 5:10]
```




    5x6 SnpArrays.SnpArray{2}:
     (true,true)   (true,true)  (false,true)  …  (true,true)   (false,true)
     (true,true)   (true,true)  (true,true)      (true,true)   (false,true)
     (false,true)  (true,true)  (true,true)      (false,true)  (true,true) 
     (false,true)  (true,true)  (true,true)      (true,true)   (false,true)
     (true,true)   (true,true)  (true,true)      (true,true)   (false,true)




```julia
# filter out rare SNPs with MAF < 0.05
hapmap[:, maf .≥ 0.05]
```




    324x12085 SnpArrays.SnpArray{2}:
     (true,true)   (false,false)  (true,true)   …  (false,true)  (false,true)
     (false,true)  (false,true)   (false,true)     (true,true)   (true,true) 
     (true,true)   (false,true)   (false,true)     (true,true)   (true,true) 
     (true,true)   (false,true)   (true,true)      (false,true)  (false,true)
     (true,true)   (false,true)   (false,true)     (true,true)   (true,true) 
     (false,true)  (true,true)    (true,true)   …  (false,true)  (false,true)
     (true,true)   (true,true)    (true,true)      (true,true)   (true,true) 
     (true,true)   (false,false)  (true,true)      (true,true)   (true,true) 
     (true,true)   (false,true)   (false,true)     (true,true)   (true,true) 
     (true,true)   (false,true)   (true,true)      (false,true)  (false,true)
     (true,true)   (false,true)   (true,true)   …  (true,true)   (true,true) 
     (true,true)   (true,true)    (false,true)     (false,true)  (false,true)
     (true,true)   (false,false)  (true,true)      (false,true)  (false,true)
     ⋮                                          ⋱                            
     (true,true)   (false,true)   (true,true)      (false,true)  (false,true)
     (true,true)   (false,false)  (true,true)      (false,true)  (false,true)
     (true,true)   (false,false)  (true,true)      (true,true)   (true,true) 
     (true,true)   (false,false)  (false,true)  …  (true,true)   (true,true) 
     (true,true)   (false,true)   (true,true)      (true,true)   (true,true) 
     (true,true)   (true,true)    (false,true)     (false,true)  (false,true)
     (true,true)   (false,true)   (false,true)     (false,true)  (false,true)
     (true,true)   (false,true)   (true,true)      (true,true)   (true,true) 
     (true,true)   (false,false)  (true,true)   …  (false,true)  (false,true)
     (true,true)   (false,true)   (false,true)     (false,true)  (false,true)
     (true,true)   (false,false)  (false,true)     (false,true)  (false,true)
     (true,true)   (false,false)  (true,true)      (true,true)   (true,true) 




```julia
# filter out individuals with genotyping success rate < 0.90
hapmap[missings_by_person / people .< 0.1, :]
```




    220x13928 SnpArrays.SnpArray{2}:
     (true,true)  (true,true)   (false,false)  …  (true,true)   (true,true)
     (true,true)  (false,true)  (false,true)      (false,true)  (true,true)
     (true,true)  (true,true)   (false,true)      (true,true)   (true,true)
     (true,true)  (true,true)   (true,true)       (true,true)   (true,true)
     (true,true)  (true,true)   (false,false)     (true,true)   (true,true)
     (true,true)  (true,true)   (false,true)   …  (true,true)   (true,true)
     (true,true)  (true,true)   (false,true)      (true,true)   (true,true)
     (true,true)  (true,true)   (false,true)      (true,true)   (true,true)
     (true,true)  (true,true)   (false,false)     (true,true)   (true,true)
     (true,true)  (true,true)   (false,true)      (true,true)   (true,true)
     (true,true)  (false,true)  (false,true)   …  (true,true)   (true,true)
     (true,true)  (true,true)   (true,true)       (true,true)   (true,true)
     (true,true)  (true,true)   (true,true)       (true,true)   (true,true)
     ⋮                                         ⋱                           
     (true,true)  (true,true)   (false,false)     (true,true)   (true,true)
     (true,true)  (true,true)   (false,false)     (true,true)   (true,true)
     (true,true)  (true,true)   (false,false)  …  (true,true)   (true,true)
     (true,true)  (true,true)   (false,true)      (true,true)   (true,true)
     (true,true)  (true,true)   (false,false)     (true,true)   (true,true)
     (true,true)  (true,true)   (false,false)     (true,true)   (true,true)
     (true,true)  (true,true)   (false,false)     (true,true)   (true,true)
     (true,true)  (true,true)   (false,true)   …  (true,true)   (true,true)
     (true,true)  (true,true)   (false,false)     (true,true)   (true,true)
     (true,true)  (true,true)   (false,true)      (true,true)   (true,true)
     (true,true)  (true,true)   (false,false)     (true,true)   (true,true)
     (true,true)  (true,true)   (false,false)     (true,true)   (true,true)



`sub()` and `slice()` create views of subarray without copying data and improve efficiency in many calculations.


```julia
mafcommon, = summarize(sub(hapmap, :, maf .≥ 0.05))
mafcommon'
```




    1x12085 Array{Float64,2}:
     0.0776398  0.324074  0.191589  …  0.310937  0.23913  0.23913  0.23913



## Assignment

It is possible to assign specific genotypes to a `SnpArray` entry.


```julia
hapmap[1, 1]
```




    (true,true)




```julia
hapmap[1, 1] = (false, true)
hapmap[1, 1]
```




    (false,true)




```julia
hapmap[1, 1] = NaN
hapmap[1, 1]
```




    (true,false)




```julia
hapmap[1, 1] = 2
hapmap[1, 1]
```




    (true,true)



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
|---|---|---|  
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




    324x13928 Array{Float64,2}:
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




    324x13928 sparse matrix with 1614876 Float32 entries:
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
countnz(isnan(hapmap)), countnz(isnan(hapmapf64))
```




    (11894,11894)



One can enforce **crude imputation** by setting the optional argument `impute=true`. Imputation is done by generating two random alleles according to the minor allele frequency. This is a neutral but not an optimal strategy, and users should impute missing genotypes by more advanced methods.


```julia
hapmapf64impute = convert(Matrix{Float64}, hapmap; impute = true)
countnz(isnan(hapmapf64impute))
```




    0



By default `convert()` translates genotypes according to the *additive* SNP model, which essentially counts the number of **minor allele** (0, 1 or 2) per genotype. Other SNP models are *dominant* and *recessive*, both in terms of the **minor allele**. When `A1` is the minor allele, genotypes are translated to real number according to

| Genotype | `SnpArray` | `model=:additive` | `model=:dominant` | `model=:recessive` |    
|---|---|---|---|---|  
| A1,A1 | (false,false) | 2 | 1 | 1 |  
| A1,A2 | (false,true) | 1 | 1 | 0 |  
| A2,A2 | (true,true) | 0 | 0 | 0 |  
| missing | (true,false) | NaN | NaN | NaN | 

When `A2` is the minor allele, genotypes are translated according to

| Genotype | `SnpArray` | `model=:additive` | `model=:dominant` | `model=:recessive` |    
|---|---|---|---|---|  
| A1,A1 | (false,false) | 0 | 0 | 0 |  
| A1,A2 | (false,true) | 1 | 1 | 0 |  
| A2,A2 | (true,true) | 2 | 1 | 1 |  
| missing | (true,false) | NaN | NaN | NaN |


```julia
[convert(Vector{Float64}, hapmap[1:10, 5]; model = :additive) convert(Vector{Float64}, hapmap[1:10, 5]; model = :dominant) convert(Vector{Float64}, hapmap[1:10, 5]; model = :recessive)]
```




    10x3 Array{Float64,2}:
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




    324x2 Array{Float64,2}:
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




    324x324 Array{Float64,2}:
     0.56716    0.0446292  0.0187812  …  0.0623061  0.0688333  0.0623358
     0.0446292  0.530456   0.0308252     0.0496909  0.0432836  0.0607344
     0.0187812  0.0308252  0.511553      0.0444371  0.0292447  0.0346823
     0.0464774  0.0352553  0.0278898     0.0575389  0.0639488  0.0577517
     0.0501686  0.0412886  0.0239638     0.0689765  0.05581    0.0631761
     0.0430517  0.0306118  0.0373675  …  0.0680426  0.0543197  0.0634806
     0.0377096  0.0212029  0.0114879     0.0426531  0.0367017  0.0355074
     0.0401055  0.0368531  0.0208137     0.0556217  0.0527508  0.0641626
     0.0289834  0.0297831  0.0159816     0.0320428  0.0440633  0.0368502
     0.0377828  0.040659   0.0248875     0.0642697  0.0552253  0.046966 
     0.0459799  0.043705   0.0219533  …  0.0564642  0.0645977  0.0587296
     0.057849   0.0373547  0.0356479     0.0640904  0.0564836  0.0696517
     0.0347474  0.042112   0.0255084     0.0551066  0.0673372  0.0608667
     ⋮                                ⋱                                 
     0.0635043  0.0527112  0.0367059     0.0830041  0.0831143  0.0933925
     0.0645941  0.0500013  0.0318383     0.0882609  0.0773871  0.0856311
     0.0656201  0.0538004  0.0371652     0.0864613  0.081387   0.0823612
     0.062374   0.0544159  0.0380025  …  0.0841856  0.0734154  0.0859109
     0.0601551  0.0539371  0.033902      0.0759062  0.0706764  0.0789132
     0.0624709  0.0580887  0.0382063     0.0848764  0.0885737  0.0796453
     0.0603521  0.0579762  0.0370625     0.0796814  0.081552   0.077682 
     0.0596907  0.0603143  0.0318682     0.0902616  0.0863425  0.0831111
     0.0640837  0.0609359  0.0311162  …  0.0938209  0.0760241  0.0830538
     0.0623061  0.0496909  0.0444371     0.606927   0.0788847  0.0867189
     0.0688333  0.0432836  0.0292447     0.0788847  0.583824   0.0792973
     0.0623358  0.0607344  0.0346823     0.0867189  0.0792973  0.576008 




```julia
# GRM using every other SNP
grm(sub(hapmap, :, 1:2:snps))
```




    324x324 Array{Float64,2}:
     0.555711   0.0411081  0.0265089  …  0.064926   0.0713428  0.0655421
     0.0411081  0.548667   0.034898      0.0554622  0.0442424  0.0542212
     0.0265089  0.034898   0.503912      0.0381145  0.0378988  0.0457407
     0.0431157  0.0444808  0.0251197     0.0491527  0.0599274  0.0541949
     0.0504286  0.0462405  0.0257578     0.0648966  0.054637   0.06293  
     0.0503159  0.0390761  0.0375003  …  0.0743448  0.0598597  0.0514303
     0.0379731  0.0254495  0.016932      0.0451145  0.0330481  0.0324201
     0.0454739  0.0364401  0.0249297     0.056584   0.0536136  0.0665555
     0.0255922  0.0218296  0.0187436     0.0319307  0.0418792  0.0342622
     0.0311334  0.0406799  0.0230293     0.0574701  0.0391224  0.0481948
     0.0470478  0.0487236  0.020016   …  0.0605443  0.062677   0.0514672
     0.060893   0.0489569  0.0379026     0.0635047  0.0620609  0.0620869
     0.0244159  0.0416684  0.0220674     0.0525605  0.0662407  0.0612487
     ⋮                                ⋱                                 
     0.0578218  0.0506541  0.0386903     0.0838834  0.0817911  0.0979834
     0.0638717  0.0581029  0.0303247     0.0861111  0.0736315  0.0840544
     0.0680155  0.0519169  0.038517      0.0816073  0.078762   0.0920692
     0.0631713  0.0561309  0.0397746  …  0.0900536  0.0740586  0.0799345
     0.060753   0.0605187  0.0400069     0.0705926  0.0671245  0.0840141
     0.0632702  0.0633616  0.0394052     0.0735999  0.0797712  0.0864385
     0.0607     0.0614583  0.030548      0.0825165  0.0894923  0.0802724
     0.0580587  0.0590514  0.0336872     0.0916351  0.0881917  0.0805812
     0.073517   0.0639522  0.0387537  …  0.0981214  0.0718733  0.0864357
     0.064926   0.0554622  0.0381145     0.599362   0.0670731  0.092225 
     0.0713428  0.0442424  0.0378988     0.0670731  0.561742   0.0721434
     0.0655421  0.0542212  0.0457407     0.092225   0.0721434  0.576505 




```julia
# MoM using all SNPs
grm(hapmap; method = :MoM)
```




    324x324 Array{Float64,2}:
     0.53945     0.0350882  0.00237734  …  0.0538645  0.0634298  0.050558 
     0.0350882   0.518194   0.0144225      0.042646   0.0395756  0.0504399
     0.00237734  0.0144225  0.499654       0.0318998  0.0205632  0.0189099
     0.0428821   0.0288295  0.022925       0.0520931  0.0695704  0.0493771
     0.045362    0.033553   0.0154853      0.0656735  0.0566986  0.0559901
     0.0318998   0.020327   0.0262315   …  0.0596509  0.0502037  0.049259 
     0.0246963   0.0113522  0.00308588     0.0314274  0.0310732  0.0244601
     0.0262315   0.0291837  0.00958082     0.038749   0.0464248  0.0472515
     0.0212717   0.0257591  0.0133597      0.0291837  0.0484324  0.0300104
     0.0210355   0.0264677  0.0105255      0.0530379  0.0503218  0.0323721
     0.0359148   0.0339073  0.00544767  …  0.0494952  0.055872   0.0502037
     0.0483143   0.0326083  0.0272943      0.0553997  0.0552816  0.0525655
     0.0307189   0.0422917  0.0245782      0.0587062  0.0703971  0.0578795
     ⋮                                  ⋱                                 
     0.0422917   0.0417012  0.0257591      0.0712237  0.0791357  0.0817337
     0.062367    0.04796    0.0270581      0.091417   0.0786634  0.0836231
     0.0544549   0.0452439  0.0250506      0.0804347  0.0812613  0.0721684
     0.0455982   0.0432364  0.024224    …  0.0727589  0.0629574  0.0709875
     0.0562263   0.0535102  0.0294199      0.0780729  0.0706332  0.0741759
     0.053274    0.0546911  0.0304827      0.0790176  0.0869296  0.0718141
     0.051857    0.0498494  0.028239       0.076892   0.082206   0.0726408
     0.0484324   0.0529198  0.0209174      0.0816156  0.0838593  0.0714599
     0.050558    0.0483143  0.0211536   …  0.0882286  0.0735855  0.0696885
     0.0538645   0.042646   0.0318998      0.561533   0.080789   0.0760654
     0.0634298   0.0395756  0.0205632      0.080789   0.566847   0.0690981
     0.050558    0.0504399  0.0189099      0.0760654  0.0690981  0.533427 



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




    (
    324x3 Array{Float64,2}:
     -38.7231  -1.2983     -7.00541  
     -32.6096  -1.21052    -3.3232   
     -23.0215  -0.505397   12.1751   
     -35.692   -2.76103    -2.40055  
     -37.1815  -0.132498   -3.66829  
     -34.9285  -1.11368     6.14167  
     -22.0323  -5.70536     2.02968  
     -30.9994  -2.28269    -0.0893283
     -22.8432  -3.76024     7.97486  
     -32.2024  -0.239253    2.91168  
     -36.344   -0.773184   -5.31525  
     -35.886   -0.807234    0.279053 
     -33.9423  -3.78982     7.35677  
       ⋮                             
     -49.1282   0.913683   10.4061   
     -46.9862  -0.9654     -0.435579 
     -48.5334  -1.05076    -0.15223  
     -49.0331   0.379279    5.65431  
     -47.8714  -0.406195   -7.14605  
     -48.2028  -1.41369    -0.564107 
     -46.7128  -3.36643    -4.44341  
     -48.9006  -1.69293     0.0467995
     -48.5574   1.34936    -1.89814  
     -50.2291   0.0865293  -1.94494  
     -48.9263  -2.06102     2.17374  
     -48.8627   0.274894    6.49518  ,
    
    13928x3 Array{Float64,2}:
      9.66817e-20   7.35949e-19   5.79015e-19
      0.00143962   -0.0042375    -0.00311816 
     -0.0183601    -0.00512036    0.00322409 
     -0.00956451   -0.004523     -0.00478078 
      0.0211999    -0.0226285     0.0110026  
     -1.82e-19     -1.35541e-18  -1.07856e-18
     -0.00230269   -0.000231224  -0.00339487 
     -0.0202126    -0.0025855     8.10915e-5 
      0.00631175   -0.0181213     0.00582407 
      0.000691273  -0.00158342   -0.0121338  
     -6.34042e-19  -3.71923e-18  -2.90818e-18
      0.0186933     7.92095e-5    0.00276918 
     -0.0127722     0.00765991    0.0134646  
      ⋮                                      
      0.000732677   0.000506129   0.00241864 
      0.000632772   0.000487763   0.00243887 
     -0.000604616  -0.000224069  -0.00294191 
      0.000769648   0.000534368   0.00250158 
      0.000410429   0.000371501   0.00266287 
     -0.00115497   -0.00172623    0.00106324 
      0.00051705    0.000728929   0.00249755 
      0.000652703   0.000748617   0.0023053  
      0.000643944  -0.000151043   0.00242307 
     -0.00149825   -0.000183435  -0.00454919 
      0.00313326   -0.00427486   -0.0152038  
     -9.09523e-5   -0.00287777    0.0037855  ,
    
    [1841.3950939952633,225.32365874997188,70.70835685208192])



To use eigen-SNPs for plotting or as covariates in GWAS, we typically scale them by their standard deviations so that they have mean zero and unit variance.


```julia
# standardize eigen-SNPs before plotting or GWAS
scale!(pcscore, 1.0 ./ √(pcvariance))
std(pcscore, 1)
```




    1x3 Array{Float64,2}:
     1.0  1.0  1.0



Internally `pca()` converts `SnpArray` to the matrix of minor allele counts. The default format is `Matrix{Float64}`, which can easily exceed memory limit. Users have several options when the default `Matrix{Float64}` cannot fit into memory.  

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
       70.7085



* Use subset of SNPs


```julia
# principal components using every other SNP capture about half the variance
srand(123)
pca(sub(hapmap, :, 1:2:snps), 3)[3]
```




    3-element Array{Float64,1}:
     926.622 
     113.188 
      36.4866



* Use sparse matrix. For large data sets with majority of rare variants, `pca_sp()` is more efficient by first converting `SnpArray` to a sparse matrix (default is `SparseMatrixCSC{Float64, Int64}`) and then computing principal components using iterative algorithms. 


```julia
# approximately same answer if we use Float16 sparse matrix
srand(123)
pca_sp(hapmap, 3, SparseMatrixCSC{Float16, UInt32})[3]
```




    3-element Array{Float64,1}:
     1841.4   
      225.31  
       70.7094




```julia
# approximately same answer if we use Int8 sparse matrix
srand(123)
pca_sp(hapmap, 3, SparseMatrixCSC{Int8, UInt32})[3]
```




    3-element Array{Float64,1}:
     1841.4   
      225.328 
       70.7119



## SnpData

`SnpData` type stores all SNP information besides genotypes. We initialize `SnpData` from Plink binary files. Note all three Plink files `.bed`, `.bim`, and `.fam` need to be present.


```julia
hapmap_snpdata = SnpData("hapmap3")
fieldnames(hapmap_snpdata)
```




    14-element Array{Symbol,1}:
     :people             
     :snps               
     :personid           
     :snpid              
     :chromosome         
     :genetic_distance   
     :basepairs          
     :allele1            
     :allele2            
     :maf                
     :minor_allele       
     :snpmatrix          
     :missings_per_person
     :missings_per_snp   




```julia
hapmap_snpdata.snpid[1:10]
```




    10-element Array{AbstractString,1}:
     "rs10458597"
     "rs12562034"
     "rs2710875" 
     "rs11260566"
     "rs1312568" 
     "rs35154105"
     "rs16824508"
     "rs2678939" 
     "rs7553178" 
     "rs13376356"




```julia
hapmap_snpdata.personid[1:10]
```




    10-element Array{AbstractString,1}:
     "NA19916"
     "NA19835"
     "NA20282"
     "NA19703"
     "NA19901"
     "NA19908"
     "NA19914"
     "NA20287"
     "NA19713"
     "NA19904"



To write `SnpData` to Plink **bed** and **bim** files, we use the `writeplink` function.


```julia
#writeplink(filename, hapmap_snpdata)
```


```julia
versioninfo()
```

    Julia Version 0.4.6
    Commit 2e358ce (2016-06-19 17:16 UTC)
    Platform Info:
      System: Darwin (x86_64-apple-darwin13.4.0)
      CPU: Intel(R) Core(TM) i7-3720QM CPU @ 2.60GHz
      WORD_SIZE: 64
      BLAS: libopenblas (USE64BITINT DYNAMIC_ARCH NO_AFFINITY Sandybridge)
      LAPACK: libopenblas64_
      LIBM: libopenlibm
      LLVM: libLLVM-3.3



```julia

```
