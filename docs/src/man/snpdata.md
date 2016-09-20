
# SnpData

`SnpData` stores all SNP information in addition to genotypes. We initialize `SnpData` from Plink binary files. Note all three Plink files `.bed`, `.bim` and `.fam` need to be present.


```julia
using SnpArrays
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



To write `SnpData` to Plink **bed** and **bim** files, we use `writeplink(filename, snpdata)`.


```julia
#writeplink(filename, hapmap_snpdata)
```
