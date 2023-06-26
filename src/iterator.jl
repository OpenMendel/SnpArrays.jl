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

function mafCount(s::SnpData,startidx=1)::Int
    n00 = 0
    n01 = 0
    n11 = 0
    for i in startidx:size(s,1)
        row = s[i,:]
        totalCount += 1
        if row[5] == 0 && row[6] == 0
            n00 += 1
        else if row[5] == 0 && row[6] == 1
            n01 += 1
        else if row[5] == 1 && row[6] == 1
            n11 += 1
        end
    end
    
    return n00,n01,n11,totalCount
end

function maf(s::SnpData, startidx=1)
    results = mafCount(s)
    if results[1] <= results[2] && results[1] <= results[3]
        return float(results[1]) / float(results[4])
    end
    if results[2] <= results[1] && results[2] <= results[1]
        return float(results[2]) / float(results[4])
    end
    if results[3] <= results[1] && results[3] <= results[2]
        return float(results[3]) / float(results[4])
    end
end