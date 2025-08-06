mutable struct SnpArrayIterator <: GeneticVariantBase.VariantIterator
   snpdata::SnpData
end

mutable struct SnpArrayIndex <: GeneticVariantBase.Variant
    index::Int
end 

@inline function Base.eltype(::Type{<:VariantIterator})
    SnpArrayIndex
end

function Base.iterate(itr::SnpArrayIterator, state=1)
    if state <= 0
        throw(BoundsError(itr, state))
    end
    if state > size(itr.snpdata.snparray,2)
        return nothing
    else
        index = SnpArrayIndex(state)
        state = state + 1
        return (index, state)
    end
end

@inline function Base.length(itr::SnpArrayIterator)
    return size(itr.snpdata.snparray, 2)
end

function GeneticVariantBase.chrom(s::SnpData, snpindex::SnpArrayIndex)::String
    result = s.snp_info[snpindex.index,:chromosome]
    return result
end

function GeneticVariantBase.pos(s::SnpData, snpindex::SnpArrayIndex)::Int
    result = s.snp_info[snpindex.index,:position]
    # println("entered pos function $result $snpindex.index")
    return result
end

function GeneticVariantBase.rsid(s::SnpData, snpindex::SnpArrayIndex)::String
    result = s.snp_info[snpindex.index,:snpid]
    return result
end

#SnpData subtype of Genetic Data

function alleles(s::SnpData, snpindex::SnpArrayIndex)::Vector{String}
    allele1 = s.snp_info[snpindex.index,:allele1]
    allele2 = s.snp_info[snpindex.index,:allele2]
    return [allele1, allele2]
end

function GeneticVariantBase.alt_allele(s::SnpData, snpindex::SnpArrayIndex)::String
    alt = s.snp_info[snpindex.index,:allele2]
    return alt
end

function GeneticVariantBase.ref_allele(s::SnpData, snpindex::SnpArrayIndex)::String
    ref = s.snp_info[snpindex.index,:allele1]
    return ref
end

struct MAFData
    maf_vector::Vector{Float64}
end

# fold into GeneticVariantBase.maf function name 

function calculate_maf_data(s::SnpData)
    maf_vector = maf(s.snparray)
    result = MAFData(maf_vector)
    return result 
end

function maf_index(maf_data::MAFData, snpindex::SnpArrayIndex)
    return maf_data.maf_vector[snpindex.index]
end


function GeneticVariantBase.maf(s::SnpData, snpindex::SnpArrayIndex)
    # maf_vector = calculate_maf_data(s)
    maf_vector = maf(s.snparray)
    return maf_vector[snpindex.index]
    # return maf_vector[snpindex.index]
end 

function GeneticVariantBase.hwepval(s::SnpData, snpindex::SnpArrayIndex)
    genotypes = s.snparray[:,snpindex.index]

    n00 = sum(genotypes .== 0x00) 
    n01 = sum(genotypes .== 0x02) 
    n11 = sum(genotypes .== 0x03) 

   pval = hwe(n00,n01,n11)
   return pval

end

 # 0 for homozygous allele 1
    # 2 Heterozygous
    # 3 homozygous allele 2
    # 1 is for missing 

function GeneticVariantBase.alt_dosages!(arr::AbstractArray{T}, s::SnpData, snpindex::SnpArrayIndex; mean_impute=true) where T <: Real
    GeneticVariantBase.alt_genotypes!(arr, s, snpindex; mean_impute=true)
    return arr 
end

# make sure you can read in all genotypes for a sample
# filtering SNPS


function GeneticVariantBase.alt_genotypes!(arr::AbstractArray{T}, s::SnpData, snpindex::SnpArrayIndex; mean_impute=true) where T <: Real
    Base.copyto!(arr, @view(s.snparray[:, snpindex.index]); impute=mean_impute, center=mean_impute, scale=mean_impute)   # change impute to mean_impute 
    return arr 
end

# are we reusing the same arr 

function GeneticVariantBase.n_samples(s::SnpData)::Int
    return size(s.snparray,1)
end 

function GeneticVariantBase.n_variants(s::SnpData)::Int
    return size(s.snparray,2)
end 