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