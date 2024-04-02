mutable struct SnpArrayIterator <: VariantIterator
   snpdata::SnpData
end

mutable struct SnpArrayIndex <: Variant
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

function iterator(s::SnpData)
    iterator = SnpArrayIterator(s)
    return iterator
end

function chrom(s::SnpData, snpindex::SnpArrayIndex)::String
    result = s.snp_info[snpindex.index,:chromosome]
    return result
end

function pos(s::SnpData, snpindex::SnpArrayIndex)::Int
    result = s.snp_info[snpindex.index,:position]
    return result
end

function rsid(s::SnpData, snpindex::SnpArrayIndex)::String
    result = s.snp_info[snpindex.index,:snpid]
    return result
end

#SnpData subtype of Genetic Data

function alleles(s::SnpData, snpindex::SnpArrayIndex)::Vector{String}
    allele1 = s.snp_info[snpindex.index,:allele1]
    allele2 = s.snp_info[snpindex.index,:allele2]
    return [allele1, allele2]
end

function alt_allele(s::SnpData, snpindex::SnpArrayIndex)::String
    alt = s.snp_info[snpindex.index,:allele2]
    return alt
end

function ref_allele(s::SnpData, snpindex::SnpArrayIndex)::String
    ref = s.snp_info[snpindex.index,:allele1]
    return ref
end

struct MAFData
    maf_vector::Vector{Float64}
end

function calculate_maf_data(s::SnpData)
    maf_vector = maf(s.snparray)
    result = MAFData(maf_vector)
    return result 
end

function maf_index(maf_data::MAFData, snpindex::SnpArrayIndex)
    return maf_data.maf_vector[snpindex.index]
end

function hwepval(s::SnpData, snpindex::SnpArrayIndex)
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

function alt_dosages!(arr::AbstractArray{T}, s::SnpData, snpindex::SnpArrayIndex) where T <: Real
    # @assert size(s.snparray) == size(arr)
    copyto!(arr, @view(s.snparray[:, snpindex.index]))
    return arr 
end

function alt_genotypes!(arr::AbstractArray{T}, s::SnpData, snpindex::SnpArrayIndex) where T <: Real
    # @assert size(s.snparray) == size(arr)
    copyto!(arr, @view(s.snparray[:, snpindex.index]))
    return arr 
end
