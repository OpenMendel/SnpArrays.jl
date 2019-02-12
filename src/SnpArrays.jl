__precompile__()

module SnpArrays

using CodecZlib, LinearAlgebra, Missings, Mmap, SparseArrays, Statistics, StatsBase
import Base: IndexStyle, convert, copyto!, eltype, getindex, setindex!, length, size
import DataFrames: DataFrame, rename!
import DelimitedFiles: readdlm, writedlm
import CSV: categorical!
import CSV # for CSV.read, to avoid clash with Base.read
import Glob: glob
import LinearAlgebra: mul!
import Statistics: mean, std, var
import StatsBase: counts
export AbstractSnpArray, SnpArray, SnpBitMatrix, SnpData
export compress_plink, decompress_plink, split_plink, merge_plink, write_plink
export counts, grm, maf, mean, minorallele, missingpos, missingrate, std, var
export ADDITIVE_MODEL, DOMINANT_MODEL, RECESSIVE_MODEL

const ADDITIVE_MODEL = Val(1)
const DOMINANT_MODEL = Val(2)
const RECESSIVE_MODEL = Val(3)

include("snparray.jl")
include("filter.jl")
include("snpdata.jl")
include("grm.jl")
include("linalg.jl")

datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

end # module
