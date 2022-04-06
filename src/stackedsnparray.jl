function StackedSnpArray(arrays::Vector{SnpArray})
    @assert length(arrays) > 0
    m = arrays[1].m
    for i in 2:length(arrays)
        @assert m == arrays[i].m
    end
    ns = Array{Int}(undef, length(arrays))
    offsets = similar(ns, length(ns) + 1)
    offset = 0
    for i in 1:length(arrays)
        ns[i] = size(arrays[i], 2)
        offsets[i] = offset
        offset += ns[i]
    end
    offsets[end] = offset
    n = offset
    StackedSnpArray(arrays, m, n, ns, offsets)
end

@inline function Base.getindex(s::StackedSnpArray, i::Integer, j::Integer)
    @boundscheck checkbounds(s, i, j)
    arrayidx = searchsortedfirst(s.offsets, j) - 1
    @inbounds inarrayidx = j - s.offsets[arrayidx]
    s.arrays[arrayidx][i, inarrayidx]
end

@inline Base.eltype(s::StackedSnpArray) = UInt8

@inline Base.length(s::StackedSnpArray) = s.m * s.n

@inline Base.size(s::StackedSnpArray) = (s.m, s.n)

@inline Base.size(s::StackedSnpArray, k::Integer) = 
    k == 1 ? s.m : 
        (k == 2 ? s.n : 
            (k > 2 ? 1 : error("Dimension k out of range")))
