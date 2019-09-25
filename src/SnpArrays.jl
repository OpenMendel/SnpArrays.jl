__precompile__()

module SnpArrays

using CodecZlib, CodecXz, CodecBzip2, CodecZstd, Distributions, TranscodingStreams
using Glob, LinearAlgebra, Missings, Mmap, SparseArrays, Statistics, StatsBase
import Base: IndexStyle, convert, copyto!, eltype, getindex, setindex!, length, size
import DataFrames: DataFrame, rename!, eachrow
import DelimitedFiles: readdlm, writedlm
import CSV: categorical!
import CSV # for CSV.read, to avoid clash with Base.read
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

include("codec.jl")
include("snparray.jl")
include("filter.jl")
include("cat.jl")
include("snpdata.jl")
include("grm.jl")
include("linalg.jl")
include("reorder.jl")

datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

end # module
