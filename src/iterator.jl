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
    if state > size(itr.snpdata.snp_info,2)
        return nothing
    else
        index = SnpArrayIndex(state)
        # don't need view just return the snparray index and the next state 
        # tuple of length 2 
        state = state + 1
        return (index, state)
    end
end

@inline function Base.length(itr::SnpArrayIterator)
    return size(itr.snpdata.snp_info, 2)
end

function iterator(s::SnpData; startidx=1)
    iterator = SnpArrayIterator(s)
end

function chrom(s::SnpData, startidx=1)::String
    result = s.snp_info[startidx,1]
    return result
end

function pos(s::SnpData, startidx=1)::Int
    result = s.snp_info[startidx,2]
    return result
end

function rsid(s::SnpData, startidx=1)::String
    result = s.snp_info[startidx,3]
    return result
end

function alleles(s::SnpData, startidx=1)::Vector{String}
    allele1 = s.snp_info[startidx,5]
    allele2 = s.snp_info[startidx,6]
    return [allele1, allele2]
end

function alt_allele(s::SnpData, startidx=1)::String
    alt = s.snp_info[startidx,6]
    return alt
end

function ref_allele(s::SnpData, startidx=1)::String
    ref = s.snp_info[startidx,5]
    return ref
end

