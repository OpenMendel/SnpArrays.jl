__precompile__()

module SnpArrays

using CodecZlib, CodecXz, CodecBzip2, CodecZstd,  TranscodingStreams
using Adapt, Glob, LinearAlgebra, LoopVectorization, Missings, Mmap, Printf
using Requires, SparseArrays, Statistics, StatsBase
import Base: IndexStyle, convert, copyto!, eltype, getindex, setindex!, length, size, wait
import DataFrames: DataFrame, rename!, eachrow
import DelimitedFiles: readdlm, writedlm
import CSV # for CSV.read, to avoid clash with Base.read
import LinearAlgebra: copytri!, mul!
import Statistics: mean, std, var
import StatsBase: counts
import SpecialFunctions: gamma_inc
import VectorizationBase: gesp
import Tables: table
export AbstractSnpArray, AbstractSnpBitMatrix, AbstractSnpLinAlg
export SnpArray, SnpBitMatrix, SnpLinAlg, SnpData
export compress_plink, decompress_plink, split_plink, merge_plink, write_plink 
export counts, grm, grm_admixture, maf, mean, minorallele, missingpos, missingrate
export std, var, vcf2plink
export counts, grm, maf, mean, minorallele, missingpos, missingrate, std, var
export vcf2plink, kinship_pruning
export ADDITIVE_MODEL, DOMINANT_MODEL, RECESSIVE_MODEL
export CuSnpArray
import VariantCallFormat: findgenokey, VCF, header

const ADDITIVE_MODEL = Val(1)
const DOMINANT_MODEL = Val(2)
const RECESSIVE_MODEL = Val(3)

include("codec.jl")
include("snparray.jl")
include("filter.jl")
include("cat.jl")
include("snpdata.jl")
include("grm.jl")
include("kinship_pruning.jl")
include("linalg_direct.jl")
include("linalg_bitmatrix.jl")
include("reorder.jl")
include("vcf2plink.jl")
include("admixture.jl")

datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

function __init__()
    @require CUDA="052768ef-5323-5732-b1bb-66c8b64840ba" include("cuda.jl")
end

end # module
