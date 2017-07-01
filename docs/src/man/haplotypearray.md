
# HaplotypeArray

The `HaplotypeArray` type is similar to [SnpArray](@ref), but with two key differences:  

* The two alleles are ordered, therefore (A1, A2) is different from (A2, A1);  
* The code `(true,false)` means the genotype (A2, A1), instead of the missing genotype in `SnpArray`.  

| Genotype | HaplotypeArray |  
|:---:|:---:|  
| A1,A1 | (false,false) |  
| A1,A2 | (false,true) |  
| A2,A1 | (true,false) |  
| A2,A2 | (true,true) |  
Each bit `true` in `HaplotypeArray` indicates a copy of the A2 allele. 

## Constructor

There are various ways to initialize a `HaplotypeArray`.  

* `HaplotypeArray` can be intialized from two `BitArray`s. Each `BitArray` indicates an A2 allele copy in the first and second positions respectively.


```julia
using SnpArrays
h = HaplotypeArray(bitrand(5, 3), bitrand(5, 3))
```




    5×3 SnpArrays.HaplotypeArray{2}:
     (true, false)  (false, true)   (true, false) 
     (true, true)   (false, false)  (true, false) 
     (true, false)  (true, true)    (true, false) 
     (true, false)  (false, false)  (false, false)
     (false, true)  (true, false)   (false, false)



* `HaplotypeArray` can be intialized from a `SnpArray`.


```julia
s = SnpArray(bitrand(5, 3), bitrand(5, 3))
h = HaplotypeArray(s)
```




    5×3 SnpArrays.HaplotypeArray{2}:
     (false, true)   (true, true)    (false, true) 
     (true, false)   (true, false)   (true, true)  
     (false, false)  (false, false)  (false, true) 
     (false, false)  (false, false)  (true, false) 
     (true, true)    (true, false)   (false, false)



This constructor does **not** copy data from `SnpArray`. Therefore both `h` and `s` points to the same patch of memory. Only the interpretation of `(true, false)` changes.


```julia
isnan(s)
```




    5×3 BitArray{2}:
     false  false  false
      true   true  false
     false  false  false
     false  false   true
     false   true  false




```julia
isnan(h)
```




    5×3 BitArray{2}:
     false  false  false
     false  false  false
     false  false  false
     false  false  false
     false  false  false



Changes to `h` also effect `s`.


```julia
h[isnan(s)] = (true, true)
isnan(s)
```




    5×3 BitArray{2}:
     false  false  false
     false  false  false
     false  false  false
     false  false  false
     false  false  false



* `HaplotypeArray(m, n)` generates an m by n `HaplotypeArray` of all A1 alleles.


```julia
HaplotypeArray(5, 3)
```




    5×3 SnpArrays.HaplotypeArray{2}:
     (false, false)  (false, false)  (false, false)
     (false, false)  (false, false)  (false, false)
     (false, false)  (false, false)  (false, false)
     (false, false)  (false, false)  (false, false)
     (false, false)  (false, false)  (false, false)



## Summary statistics

`summarize` when applied to a `HaplotypeArray` only computes  

* `maf`: minor allele frequencies, taking into account of missingness.  
* `minor_allele`: a `BitVector` indicating the minor allele for each SNP.   `minor_allele[j]==true` means A1 is the minor allele for SNP j; `minor_allele[j]==false` means A2 is the minor allele for SNP j.  


```julia
maf, minor_allele = summarize(h)
```




    ([0.5, 0.4, 0.4], Bool[false, true, true])



## Subsetting and assignment

Subsetting and assignment work the same as [SnpArray](@ref).

## Copy and convert

Copying or converting a `HaplotypeArray` or slices of it to numeric arrays of **minor allele counts** is similar to `SnpArray` with the exception there are no missing genotypes in `HaplotypeArray`. So the keyword `impute` is not relevant anymore.


```julia
# convert to Matrix{Float64}
h_f64 = convert(Matrix{Float64}, h)
```




    5×3 Array{Float64,2}:
     1.0  0.0  1.0
     2.0  0.0  0.0
     0.0  2.0  1.0
     0.0  2.0  0.0
     2.0  0.0  2.0



By default `convert` translates genotypes according to the *additive* SNP model, which essentially counts the number of **minor allele** (0, 1 or 2) per genotype. Other SNP models are *dominant* and *recessive*, both in terms of the **minor allele**. When `A1` is the minor allele, genotypes are translated to real number according to

| Genotype | `HaplotypeArray` | `model=:additive` | `model=:dominant` | `model=:recessive` |    
|:---:|:---:|:---:|:---:|:---:|  
| A1,A1 | 00 | 2 | 1 | 1 |  
| A1,A2 | 01 | 1 | 1 | 0 |  
| A2,A1 | 10 | 1 | 1 | 0 |  
| A2,A2 | 11 | 0 | 0 | 0 |  

When `A2` is the minor allele, genotypes are translated according to

| Genotype | `HaplotypeArray` | `model=:additive` | `model=:dominant` | `model=:recessive` |    
|:---:|:---:|:---:|:---:|:---:|  
| A1,A1 | 00 | 0 | 0 | 0 |  
| A1,A2 | 01 | 1 | 1 | 0 |  
| A2,A1 | 01 | 1 | 1 | 0 |  
| A2,A2 | 11 | 2 | 1 | 1 |  


```julia
[convert(Vector{Float64}, h[1:5, 3]; model = :additive) convert(Vector{Float64}, h[1:5, 3]; model = :dominant) convert(Vector{Float64}, h[1:5, 3]; model = :recessive)]
```




    5×3 Array{Float64,2}:
     1.0  1.0  0.0
     0.0  0.0  0.0
     1.0  1.0  0.0
     0.0  0.0  0.0
     2.0  1.0  1.0



By default `convert` does **not** center and scale genotypes. Setting the optional arguments `center=true, scale=true` centers genotypes at 2MAF and scales them by $[2 \cdot \text{MAF} \cdot (1 - \text{MAF})]^{-1/2}$. Mono-allelic SNPs (MAF=0) are not scaled.


```julia
[convert(Vector{Float64}, h[:, 3]) convert(Vector{Float64}, h[:, 3]; center = true, scale = true)]
```




    5×2 Array{Float64,2}:
     1.0   0.288675
     0.0  -1.1547  
     1.0   0.288675
     0.0  -1.1547  
     2.0   1.73205 



`copy!` is the in-place version of `convert()`. Options such as GWAS loop over SNPs and perform statistical anlaysis for each SNP. This can be achieved by


```julia
g = zeros(size(h, 1))
for j = 1:size(h, 2)
    copy!(g, h[:, j]; model = :additive)
    # do statistical anlaysis
end
```
