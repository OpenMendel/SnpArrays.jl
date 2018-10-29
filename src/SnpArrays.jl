__precompile__()

module SnpArrays
    using LinearAlgebra, Missings, Mmap, OffsetArrays, SparseArrays, Statistics, StatsBase
    import Base: IndexStyle, convert, copyto!, eltype, getindex, setindex!, length, size
    import Statistics: mean, std, var
    import StatsBase: counts
    export SnpArray
    export bedvals, counts, grm, maf, mean, minorallele, missingpos, missingrate, std, var
    
    include("snparray.jl")

"""
    BEDvals

`Vector{Union{UInt8,Missings.Missing}}` of the possible values in a BEDFile
"""
    const bedvals = OffsetArray{Union{Int8,Missing}}(undef, 0:3)
    bedvals[0] = 0
    bedvals[2] = 1
    bedvals[3] = 2

    datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

end # module
