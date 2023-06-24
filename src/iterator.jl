mutable struct SnpArrayIterator <: VariantIterator
   snpdata::SnpData
end

struct SnpArrayIndex <: Variant
    index::Int
end 

@inline function Base.eltype(::Type{<:VariantIterator})
    SnpArray
end

function Base.iterate(itr::SnpArrayIterator, state=itr.index)
    if state > size(itr.snpdata.snp_info,2)
        return nothing
    else
        state = state + 1
        result = (view(itr.snpdata.snp_info, :, state), state+1)
        return result
    end
end

@inline function Base.length(itr::SnpArrayIterator)
    return size(itr.snpdata.snp_info, 2)
end

@inline function snparrayiterator(snpdata::SnpData)
    SnpArrayIterator(snpdata)
end

function iterator(snpdata::SnpData; startidx=1)::SnpArrayIterator
    if startidx > size(snpdata,2)
        return nothing
    else
        SnpArrayIndex(startidx)
        return SnpArrayIterator(snpdata)
    end
end

function iterator(s::SnpData; startidx=1)
    iterator = SnpArrayIterator(s)
end

function chrom(s::SnpData, startidx=1)::String
    result = view(s.snp_info,startidx,1)
    return result
end

function pos(s::SnpData, startidx=1)::Int
    result = view(s.snp_info,startidx,2)
    return result
end

function rsid(s::SnpData, startidx=1)::String
    result = view(s.snp_info,startidx,3)
    return result
end

function alleles(s::SnpData, startidx=1)::Vector{String}
    allele1 = view(s.snp_info,startidx,5)
    allele2 = view(s.snp_info,startidx,6)
    return [allele1, allele2]
end

function alt_allele(s::SnpData, startidx=1)::String
    alt = view(s.snp_info,startidx,6)
    return alt
end

function ref_allele(s::SnpData, startidx=1)::String
    ref = view(s.snp_info,v,5)
    return ref
end