__precompile__()

module SnpArrays

using LinearAlgebra, Missings, Mmap, OffsetArrays, SparseArrays, Statistics, StatsBase
import Base: IndexStyle, convert, copyto!, eltype, getindex, setindex!, length, size
import Statistics: mean, std, var
import StatsBase: counts
export AbstractSnpArray, SnpArray
export counts, grm, maf, mean, minorallele, missingpos, missingrate, std, var
export ADDITIVE_MODEL, DOMINANT_MODEL, RECESSIVE_MODEL

const ADDITIVE_MODEL = Val(1)
const DOMINANT_MODEL = Val(2)
const RECESSIVE_MODEL = Val(3)

include("snparray.jl")

datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

end # module
