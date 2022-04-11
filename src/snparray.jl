"""
    SnpArray
    SnpArray(bednm, m)
    SnpArray(plknm)
    SnpArray(undef, m, n)

Raw .bed file as a shared, memory-mapped Matrix{UInt8}.  The number of rows, `m`
is stored separately because it is not uniquely determined by the size of the `data` field.
"""
struct SnpArray <: AbstractMatrix{UInt8}
    data::Matrix{UInt8}
    columncounts::Matrix{Int}
    rowcounts::Matrix{Int}
    m::Int
end

"""
    StackedSnpArray(s::Vector{SnpArray})

Stacked SnpArray for unified indexing
"""
struct StackedSnpArray <: AbstractMatrix{UInt8} # details in stackedsnparray.jl
    arrays::Vector{SnpArray}
    m::Int
    n::Int
    ns::Vector{Int}
    offsets::Vector{Int} # 0-based
end

const AbstractSnpArray = Union{SnpArray, SubArray{UInt8, 1, SnpArray}, SubArray{UInt8, 2, SnpArray}, 
    StackedSnpArray, SubArray{UInt8, 1, StackedSnpArray}, SubArray{UInt8, 2, StackedSnpArray}}

function SnpArray(bednm::AbstractString, m::Integer, args...; kwargs...)
    checkplinkfilename(bednm, "bed")
    data = makestream(bednm, args...; kwargs...) do io
        read(io, UInt16) == 0x1b6c || throw(ArgumentError("wrong magic number in file $bednm"))
        read(io, UInt8) == 0x01 || throw(ArgumentError(".bed file, $bednm, is not in correct orientation"))
        if endswith(bednm, ".bed")
            return Mmap.mmap(io)
        else
            return read(io)
        end
    end
    drows = (m + 3) >> 2   # the number of rows in the Matrix{UInt8}
    n, r = divrem(length(data), drows)
    iszero(r) || throw(ArgumentError("filesize of $bednm is not a multiple of $drows"))
    SnpArray(reshape(data, (drows, n)), zeros(Int, (4, n)), zeros(Int, (4, m)), m)
end

function SnpArray(
    bednm::AbstractString, 
    args...; 
    famnm::Union{AbstractString, Nothing}=nothing, 
    kwargs...)
    checkplinkfilename(bednm, "bed")
    # check user supplied fam filename
    if famnm === nothing
        famnm = replace(bednm, ".bed" => ".fam") 
        isfile(famnm) || throw(ArgumentError("fam file not found"))
    else # user supplied fam file name
        checkplinkfilename(famnm, "fam")
        isfile(famnm) || throw(ArgumentError("file $famnm not found"))
    end
    m = makestream(famnm) do stream
        countlines(stream)
    end
    SnpArray(bednm, m, args...; kwargs...)
end

function SnpArray(::UndefInitializer, m::Integer, n::Integer)
    SnpArray(Matrix{UInt8}(undef, ((m + 3) >> 2, n)), zeros(Int, (4, n)), zeros(Int, (4, m)), m)
end

function SnpArray(file::AbstractString, s::SnpArray)
    makestream(file, "w+") do io
        write(io, 0x1b6c)
        write(io, 0x01)
        write(io, s.data)
    end
    SnpArray(file, s.m, "r+")
end

function SnpArray(file::AbstractString, m::Integer, n::Integer)
    makestream(file, "w+") do io
        write(io, 0x1b6c)
        write(io, 0x01)
        write(io, fill(0x00, ((m + 3) >> 2, n)))
    end
    SnpArray(file, m, "r+")
end

StatsBase.counts(s::AbstractSnpArray; dims=:) = _counts(s, dims)

function _counts(s::AbstractSnpArray, dims::Integer)
    if isone(dims)
        cc = (typeof(s) == SnpArray) ? s.columncounts : zeros(Int, (4, size(s, 2)))
        if all(iszero, cc)
            m, n = size(s)
            @inbounds for j in 1:n
                for i in 1:m
                    cc[s[i, j] + 1, j] += 1            
                end
            end
        end
        return cc
    elseif dims == 2
        rc = (typeof(s) == SnpArray) ? s.rowcounts : zeros(Int, (4, size(s, 1)))
        if all(iszero, rc)
            m, n = size(s)
            @inbounds for j in 1:n
                for i in 1:m
                    rc[s[i, j] + 1, i] += 1
                end
            end
        end
        return rc
    else
        throw(ArgumentError("counts(s::SnpArray, dims=k) only defined for k = 1 or 2"))
    end
end        

_counts(s::AbstractSnpArray, ::Colon) = sum(_counts(s, 1), dims=2)

function Base.getindex(s::SnpArray, i::Int)  # Linear indexing
    d, r = divrem(i - 1, s.m)
    s[r + 1, d + 1]
end

@inline function Base.getindex(s::SnpArray, i::Integer, j::Integer)
    @boundscheck checkbounds(s, i, j)
    ip3 = i + 3
    (s.data[ip3 >> 2, j] >> ((ip3 & 0x03) << 1)) & 0x03
end

function Base.setindex!(s::SnpArray, x::UInt8, i::Int)  # Linear indexing
    d, r = divrem(i - 1, s.m)
    Base.setindex!(s, x, r + 1, d + 1)
end

@inline function Base.setindex!(s::SnpArray, x::UInt8, i::Integer, j::Integer)
    @boundscheck checkbounds(s, i, j)
    ip3 = i + 3
    shft = (ip3 & 0x03) << 1
    mask = ~(0x03 << shft)
    s.data[ip3 >> 2, j] = (s.data[ip3 >> 2, j] & mask) | (x << shft)
    x
end

Base.eltype(s::SnpArray) = UInt8

Base.length(s::SnpArray) = s.m * size(s.data, 2)

Statistics.mean(s::AbstractSnpArray; dims=:, model=ADDITIVE_MODEL) = _mean(s, dims, model)

function _mean(s::AbstractSnpArray,  dims::Integer, ::typeof(ADDITIVE_MODEL))
    m, n = size(s)
    if isone(dims)
        cc = _counts(s, 1)   # need to use extractor to force evaluation if needed
        means = Matrix{Float64}(undef, (1, n))
        @inbounds for j in 1:n
            means[j] = (cc[3, j] + 2cc[4, j]) / (cc[1, j] + cc[3, j] + cc[4, j])
        end
        return means
    elseif dims == 2
        rc = _counts(s, 2)
        means = Matrix{Float64}(undef, (m, 1))
        @inbounds for i in 1:m
            means[i] = (rc[3, i] + 2rc[4, i]) / (rc[1, i] + rc[3, i] + rc[4, i])
        end
        return means
    else
        throw(ArgumentError("mean(s::SnpArray, dims=k) only defined for k = 1 or 2"))
    end
end

function _mean(s::AbstractSnpArray,  dims::Integer, ::typeof(DOMINANT_MODEL))
    m, n = size(s)
    if isone(dims)
        cc = _counts(s, 1)   # need to use extractor to force evaluation if needed
        means = Matrix{Float64}(undef, (1, n))
        @inbounds for j in 1:n
            means[j] = (cc[3, j] + cc[4, j]) / (cc[1, j] + cc[3, j] + cc[4, j])
        end
        return means
    elseif dims == 2
        rc = _counts(s, 2)
        means = Matrix{Float64}(undef, (m, 1))
        @inbounds for i in 1:m
            means[i] = (rc[3, i] + rc[4, i]) / (rc[1, i] + rc[3, i] + rc[4, i])
        end
        return means
    else
        throw(ArgumentError("mean(s::SnpArray, dims=k) only defined for k = 1 or 2"))
    end
end

function _mean(s::AbstractSnpArray,  dims::Integer, ::typeof(RECESSIVE_MODEL))
    m, n = size(s)
    if isone(dims)
        cc = _counts(s, 1)   # need to use extractor to force evaluation if needed
        means = Matrix{Float64}(undef, (1, n))
        @inbounds for j in 1:n
            means[j] = cc[4, j] / (cc[1, j] + cc[3, j] + cc[4, j])
        end
        return means
    elseif dims == 2
        rc = _counts(s, 2)
        means = Matrix{Float64}(undef, (m, 1))
        @inbounds for i in 1:m
            means[i] = rc[4, i] / (rc[1, i] + rc[3, i] + rc[4, i])
        end
        return means
    else
        throw(ArgumentError("mean(s::SnpArray, dims=k) only defined for k = 1 or 2"))
    end
end

function _mean(s::AbstractSnpArray, ::Colon, ::typeof(ADDITIVE_MODEL))
    rc = _counts(s, 2)
    @views (sum(rc[3, :]) + 2sum(rc[4, :])) /  sum(rc[[1, 3, 4], :])
end

function _mean(s::AbstractSnpArray, ::Colon, ::typeof(DOMINANT_MODEL))
    rc = _counts(s, 2)
    @views (sum(rc[3, :]) + sum(rc[4, :])) /  sum(rc[[1, 3, 4], :])
end

function _mean(s::AbstractSnpArray, ::Colon, ::typeof(RECESSIVE_MODEL))
    rc = _counts(s, 2)
    @views sum(rc[4, :]) /  sum(rc[[1, 3, 4], :])
end

Base.size(s::SnpArray) = s.m, size(s.data, 2)

Base.size(s::SnpArray, k::Integer) = 
k == 1 ? s.m : k == 2 ? size(s.data, 2) : k > 2 ? 1 : error("Dimension k out of range")

# TODO: need to implement `var` for different SNP models

Statistics.var(s::AbstractSnpArray; corrected::Bool=true, mean=nothing, dims=:) = _var(s, corrected, mean, dims)

function _var(s::AbstractSnpArray, corrected::Bool, mean, dims::Integer)
    m, n = size(s)
    means = something(mean, Statistics.mean(s, dims=dims))
    if isone(dims)
        cc = _counts(s, 1)
        vars = Matrix{Float64}(undef, (1, n))
        for j in 1:n
            mnj = means[j]
            vars[j] = (abs2(mnj) * cc[1, j] + abs2(1.0 - mnj) * cc[3, j] + abs2(2.0 - mnj) * cc[4,j]) /
            (cc[1, j] + cc[3, j] + cc[4, j] - (corrected ? 1 : 0))
        end
        return vars
    elseif dims == 2
        rc = _counts(s, 2)
        vars = Matrix{Float64}(undef, (m, 1))
        for i in 1:m
            mni = means[i]
            vars[i] = (abs2(mni) * rc[1, i] + abs2(1.0 - mni) * rc[3, i] + abs2(2.0 - mni) * rc[4, i]) /
            (rc[1, i] + rc[3, i] + rc[4, i] - (corrected ? 1 : 0))
        end
        return vars
    end
    throw(ArgumentError("var(s::SnpArray, dims=k) only defined for k = 1 or 2"))
end

"""
    maf!(out, s)

Populate `out` with minor allele frequencies of SnpArray `s`.
"""
function maf!(out::AbstractVector{T}, s::AbstractSnpArray) where T <: AbstractFloat
    cc = _counts(s, 1)
    @inbounds for j in 1:size(s, 2)
        out[j] = (cc[3, j] + 2cc[4, j]) / 2(cc[1, j] + cc[3, j] + cc[4, j])
        (out[j] > 0.5) && (out[j] = 1 - out[j])
    end
    out
end
maf(s::AbstractSnpArray) = maf!(Vector{Float64}(undef, size(s, 2)), s)

"""
    minorallele(out, s)

Populate `out` with minor allele indicators. `out[j] == true` means A2 is the minor 
allele of `j`th column; `out[j] == false` means A1 is the minor allele.
"""
function minorallele!(out::AbstractVector{Bool}, s::AbstractSnpArray)
    cc = _counts(s, 1)
    @inbounds for j in 1:size(s, 2)
        out[j] = cc[1, j] > cc[4, j]
    end
    out
end

"""
    minorallele(s)

Calculate minor allele indicators. `out[j] == true` means A2 is the minor 
allele of `j`th column; `out[j] == false` means A1 is the minor allele.
"""
minorallele(s::AbstractSnpArray) = minorallele!(BitVector(undef, size(s, 2)), s)

@inline function convert(::Type{T}, x::UInt8, ::Val{1}) where T <: AbstractFloat
    iszero(x) ? zero(T) : isone(x) ? T(NaN) : T(x - 1)
end

@inline function convert(::Type{T}, x::UInt8, ::Val{2}) where T <: AbstractFloat
    iszero(x) ? zero(T) : isone(x) ? T(NaN) : one(T)
end

@inline function convert(::Type{T}, x::UInt8, ::Val{3}) where T <: AbstractFloat
    (iszero(x) || x == 2) ? zero(T) : isone(x) ? T(NaN) : one(T)
end

"""
    Base.copyto!(v, s, model=ADDITIVE_MODEL, center=false, scale=false, impute=false)

Copy SnpArray `s` to numeric vector or matrix `v`.

# Arguments
- `model::Union{Val{1}, Val{2}, Val{3}}=ADDITIVE_MODEL`: `ADDITIVE_MODEL` (default), `DOMINANT_MODEL`, or `RECESSIVE_MODEL`.  
- `center::Bool=false`: center column by mean.
- `scale::Bool=false`: scale column by theoretical variance.
- `impute::Bool=false`: impute missing values by column mean.
"""
function Base.copyto!(
    v::AbstractVecOrMat{T}, 
    s::AbstractSnpArray;
    model::Union{Val{1}, Val{2}, Val{3}} = ADDITIVE_MODEL,
    center::Bool = false,
    scale::Bool = false,
    impute::Bool = false
    ) where T <: AbstractFloat
    m, n = size(s, 1), size(s, 2)
    # no center, scale, or impute
    if !center && !scale && !impute
        @inbounds for j in 1:n
            @simd for i in 1:m
                v[i, j] = SnpArrays.convert(T, s[i, j], model)
            end
        end
        return v
    end
    # center, scale, impute
    @inbounds for j in 1:n
        μj, mj = zero(T), 0
        @simd for i in 1:m
            vij = SnpArrays.convert(T, s[i, j], model)
            v[i, j] = vij
            μj += isnan(vij) ? zero(T) : vij
            mj += isnan(vij) ? 0 : 1
        end
        μj /= mj
        σj = model == ADDITIVE_MODEL ? sqrt(μj * (1 - μj / 2)) : sqrt(μj * (1 - μj))
        @simd for i in 1:m
            impute && isnan(v[i, j]) && (v[i, j] = μj)
            center && (v[i, j] -= μj)
            scale && σj > 0 && (v[i, j] /= σj)
        end
    end
    return v
end


"""
    Base.convert(t, s, model=ADDITIVE_MODEL, center=false, scale=false, impute=false)

Convert a SnpArray `s` to a numeric vector or matrix of same shape as `s`.

# Arguments
- `t::Type{AbstractVecOrMat{T}}`: Vector or matrix type.
- `model::Union{Val{1}, Val{2}, Val{3}}=ADDITIVE_MODEL`: `ADDITIVE_MODEL` (default), `DOMINANT_MODEL`, or `RECESSIVE_MODEL`.  
- `center::Bool=false`: center column by mean.
- `scale::Bool=false`: scale column by theoretical variance.
- `impute::Bool=false`: impute missing values by column mean.
"""
function Base.convert(
    ::Type{T},
    s::AbstractSnpArray;
    kwargs...) where T <: Array
    T(s; kwargs...)
end
Array{T,N}(s::AbstractSnpArray; kwargs...) where {T,N} = 
    copyto!(Array{T,N}(undef, size(s)), s; kwargs...)


"""
    missingpos(s::SnpArray)

Return a `SparseMatrixCSC{Bool,Int32}` of the same size as `s` indicating the positions with missing data.
"""
function missingpos(s::AbstractSnpArray)
    m, n = size(s)
    colptr = sizehint!(Int32[1], n + 1)
    rowval = Int32[]
    @inbounds for j in 1:n
        msngpos = Int32[]
        for i in 1:m
            isone(s[i, j]) && push!(msngpos, i)
        end
        append!(rowval, msngpos)
        push!(colptr, colptr[end] + length(msngpos))
    end
    SparseMatrixCSC(m, n, colptr, rowval, fill(true, length(rowval)))
end

function _missingrate!(out::AbstractVector{<:AbstractFloat}, s::AbstractSnpArray,  dims::Integer)
    m, n = size(s)
    if isone(dims)
        cc = _counts(s, 1)   # need to use extractor to force evaluation if needed
        @inbounds for j in 1:n
            out[j] = cc[2, j] / m
        end
    elseif dims == 2
        rc = _counts(s, 2)
        @inbounds for i in 1:m
            out[i] = rc[2, i] / n
        end
    else
        throw(ArgumentError("_missingrate(out, s::SnpArray, dims=k) only defined for k = 1 or 2"))
    end
    out
end

function missingrate(s::AbstractSnpArray, dims::Integer)
    if isone(dims)
        return _missingrate!(Vector{Float64}(undef, size(s, 2)), s, 1)
    elseif dims == 2 
        return _missingrate!(Vector{Float64}(undef, s.m), s, 2)
    else
        throw(ArgumentError("missingrate(s::SnpArray, dims=k) only defined for k = 1 or 2"))
    end
end
