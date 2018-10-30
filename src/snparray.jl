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

AbstractSnpArray = Union{SnpArray, SubArray{UInt8, 1, SnpArray}, SubArray{UInt8, 2, SnpArray}}

function SnpArray(bednm::AbstractString, m::Integer, args...; kwargs...)
    data = open(bednm, args...; kwargs...) do io
        read(io, UInt16) == 0x1b6c || throw(ArgumentError("wrong magic number in file $bednm"))
        read(io, UInt8) == 0x01 || throw(ArgumentError(".bed file, $bednm, is not in correct orientation"))
        Mmap.mmap(io)
    end
    drows = (m + 3) >> 2   # the number of rows in the Matrix{UInt8}
    n, r = divrem(length(data), drows)
    iszero(r) || throw(ArgumentError("filesize of $bednm is not a multiple of $drows"))
    SnpArray(reshape(data, (drows, n)), zeros(Int, (4, n)), zeros(Int, (4, m)), m)
end

SnpArray(nm::AbstractString, args...; kwargs...) = 
    SnpArray(nm, countlines(string(splitext(nm)[1], ".fam")), args...; kwargs...)

function SnpArray(::UndefInitializer, m::Integer, n::Integer)
    SnpArray(Matrix{UInt8}(undef, ((m + 3) >> 2, n)), zeros(Int, (4, n)), zeros(Int, (4, m)), m)
end

function SnpArray(file::AbstractString, s::SnpArray)
    open(file, "w+") do io
        write(io, 0x1b6c)
        write(io, 0x01)
        write(io, s.data)
    end
    SnpArray(file, s.m, "r+")
end

function SnpArray(file::AbstractString, m::Integer, n::Integer)
    open(file, "w+") do io
        write(io, 0x1b6c)
        write(io, 0x01)
        write(io, fill(0x00, ((m + 3) >> 2, n)))
    end
    SnpArray(file, m, "r+")
end

StatsBase.counts(s::SnpArray; dims=:) = _counts(s, dims)

function _counts(s::SnpArray, dims::Integer)
    if isone(dims)
        cc = s.columncounts
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
        rc = s.rowcounts
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

_counts(s::SnpArray, ::Colon) = sum(_counts(s, 1), dims=2)

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

Statistics.mean(s::SnpArray; dims=:) = _mean(s, dims)

function _mean(s::SnpArray,  dims::Integer)
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

function _mean(s::SnpArray, ::Colon)
    rc = _counts(s, 2)
    (sum(view(rc, 3, :)) + 2sum(view(rc, 4, :))) / sum(view(rc, [1, 3, 4], :))
end

Base.size(s::SnpArray) = s.m, size(s.data, 2)

Base.size(s::SnpArray, k::Integer) = 
k == 1 ? s.m : k == 2 ? size(s.data, 2) : k > 2 ? 1 : error("Dimension k out of range")

Statistics.var(s::SnpArray; corrected::Bool=true, mean=nothing, dims=:) = _var(s, corrected, mean, dims)

function _var(s::SnpArray, corrected::Bool, mean, dims::Integer)
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
function maf!(out::AbstractVector{T}, s::SnpArray) where T <: AbstractFloat
    cc = _counts(s, 1)
    @inbounds for j in 1:size(s, 2)
        out[j] = (cc[3, j] + 2cc[4, j]) / 2(cc[1, j] + cc[3, j] + cc[4, j])
        (out[j] > 0.5) && (out[j] = 1 - out[j])
    end
    out
end
maf(s::SnpArray) = maf!(Vector{Float64}(undef, size(s, 2)), s)

"""
    minorallele(out, s)

Populate `out` with minor allele indicators. `out[j] == true` means A2 is the minor 
allele of `j`th column; `out[j[ == false` means A1 is the minor allele.
"""
function minorallele!(out::AbstractVector{Bool}, s::SnpArray)
    cc = _counts(s, 1)
    @inbounds for j in 1:size(s, 2)
        out[j] = cc[1, j] > cc[4, j]
    end
    out
end
minorallele(s::SnpArray) = minorallele!(Vector{Bool}(undef, size(s, 2)), s)

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
- `impute::Bool=falase`: impute missing values by column mean.

# Example

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
- `impute::Bool=falase`: impute missing values by column mean.
"""
function Base.convert(
    ::Type{T},
    s::AbstractSnpArray;
    kwargs...) where T <: Array
    T(s; kwargs...)
end
Array{T,N}(s::AbstractSnpArray; kwargs...) where {T,N} = copyto!(Array{T,N}(undef, size(s)), s; kwargs...)


"""    
    outer(s::SnpArray, colinds)
    outer(s::SnpArray)

Return the "outer product", `f * f'` using the `Float32[0, NaN, 1, 2]` encoding of `s`.

The `colinds` argument, when given, causes the operation to be performed on that subset
of the columns.
"""
function outer(s::SnpArray, colinds::AbstractVector{<:Integer})
    m = size(s, 1)
    outer!(Symmetric(zeros(Float32, (m, m))), s, colinds)
end    
outer(s::SnpArray) = outer(s, 1:size(s, 2))

"""
    outer!(sy::Symmetric, s::SnpArray, colinds)

Update `sy` with the sum of the outer products of the columns in `colind` from `f`    
"""
function outer!(sy::Symmetric{T}, s::SnpArray, colinds::AbstractVector{<:Integer}) where T
    tempv = Vector{T}(undef, s.m)
    for j in colinds
        LinearAlgebra.BLAS.syr!(sy.uplo, one(T), copyto!(tempv, s, j), sy.data)
    end
    sy
end

"""
    missingpos(s::SnpArray)

Return a `SparseMatrixCSC{Bool,Int32}` of the same size as `s` indicating the positions with missing data.
"""
function missingpos(s::SnpArray)
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

function _missingrate!(out::AbstractVector{<:AbstractFloat}, s::SnpArray,  dims::Integer)
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

function missingrate(s::SnpArray, dims::Integer)
    if isone(dims)
        return _missingrate!(Vector{Float64}(undef, size(s, 2)), s, 1)
    elseif dims == 2 
        return _missingrate!(Vector{Float64}(undef, s.m), s, 2)
    else
        throw(ArgumentError("missingrate(s::SnpArray, dims=k) only defined for k = 1 or 2"))
    end
end

"""
    grm(s, method=:GRM, minmaf=0.01, colinds=nothing)

Compute empirical kinship matrix from a SnpArray `s`. Missing genotypes are imputed
on the fly by column mean.

# Arguments
- `method::Symbol=:GRM`: `:GRM` (default), `:MoM`, or `Robust`.
- `minmaf::Real=0.01`: columns with MAF `<minmaf` are excluded; default 0.01.
- `cinds=nothing`: indices or mask of columns to be used for calculating GRM.
- `t::Type{T}=Float64`: Float type for calculating GRM; default is `Float64`.
"""
function grm(
    s::SnpArray;
    method::Symbol = :GRM,
    minmaf::Real = 0.01,
    cinds::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    t::Type{T} = Float64
    ) where T <: AbstractFloat
    mf = maf(s)
    colinds = something(cinds, mf .≥ minmaf)
    @views grm(s[:, colinds], mf[colinds], method, t)
end # function grm

function grm(
    s::AbstractSnpArray,
    maf::AbstractVector,
    method::Symbol = :GRM,
    ::Type{T} = Float64
    ) where T <: AbstractFloat
    m, n = size(s)
    G = Mmap.mmap(Matrix{T}, m, n)
    if method == :GRM
        Base.copyto!(G, s, model=ADDITIVE_MODEL, impute=true, center=true, scale=true)
        Φ = G * transpose(G)
        Φ ./= 2n
    elseif method == :MoM
        Base.copyto!(G, s, model=ADDITIVE_MODEL, impute=true)
        G .-= 1
        Φ = G * transpose(G)
        c = sum(x -> abs2(x) + abs2(1 - x), maf)
        shft, scal = n / 2 - c, 1 / (n - c)
        @inbounds @simd for i in eachindex(Φ)
            Φ[i] = (Φ[i] / 2 + shft) * scal
        end
    elseif method == :Robust
        Base.copyto!(G, s, model=ADDITIVE_MODEL, center=true, impute=true)
        scal = sum(x -> 4x * (1 - x), maf)
        Φ = G * transpose(G)
        Φ ./= scal
    else
        throw(ArgumentError("method should be :GRM, :MoM, or :Robust"))
    end
    Φ
end # function grm

"""
    SnpArrays.filter(s[, min_success_rate_per_row, min_success_rate_per_col, maxiters])

Filter a SnpArray by genotyping success rate.

# Input
- `A`: a SnpArray.
- `min_success_rate_per_row`: threshold for SNP genotyping success rate.
- `min_success_rate_per_col`: threshold for person genotyping success rate.
- `maxiters`: maximum number of filtering iterations.

# Output
- `rmask`: BitVector indicating remaining rows.
- `cmask`: BitVector indicating remaining cols.
"""
function filter(
    s::SnpArray, 
    min_success_rate_per_row::Real = 0.98,
    min_success_rate_per_col::Real = 0.98,
    maxiters::Integer = 5)
    m, n = size(s)
    rc, cc = zeros(Int, m), zeros(Int, n)
    rmask, cmask = trues(m), trues(n)
    rmiss, cmiss = 1 - min_success_rate_per_row, 1 - min_success_rate_per_col
    for iter in 1:maxiters
        fill!(rc, 0)
        fill!(cc, 0)
        @inbounds for j in 1:n
            cmask[j] || continue
            for i in 1:m
                rmask[i] && s[i, j] == 0x01 && (rc[i] += 1; cc[j] += 1)
            end
        end
        rows, cols = count(rmask), count(cmask)
        @inbounds for j in 1:n
            cmask[j] = cmask[j] && cc[j] < cmiss * rows
        end
        @inbounds for i in 1:m
            rmask[i] = rmask[i] && rc[i] < rmiss * cols
        end
        count(cmask) == cols && count(rmask) == rows && break
        iter == maxiters && (@warn "success rate not satisfied; consider increase maxiters")
    end
    rmask, cmask
end

"""
    SnpArrays.filter(src, rowinds, colinds; des = src * ".filtered")

Filter `src` Plink files according to row indices `rowinds` and column indices 
`colinds` and write to a new set of Plink files `des`.

# Input
- `src`: source Plink file name without suffix ".bed", ".fam" or ".bim".
- `rowinds`: row indices.
- `colinds`: column indices.

# Keyword arguments
- `des`: output Plink file name; defualt is `src * ".filtered"`.
"""
function filter(
    src::AbstractString, 
    rowinds::AbstractVector{<:Integer},
    colinds::AbstractVector{<:Integer};
    des::AbstractString = src * ".filtered")
    # check source plink files
    isfile(src * ".bed") || throw(ArgumentError("$src.bed file not found"))
    isfile(src * ".bim") || throw(ArgumentError("$src.bim file not found"))
    isfile(src * ".fam") || throw(ArgumentError("$src.fam file not found"))
    # create row and column masks
    if eltype(rowinds) == Bool
        rmask = rowinds
    else
        rmask = falses(countlines(src * ".fam"))
        rmask[rowinds] .= true
    end
    if eltype(colinds) == Bool
        cmask = colinds
    else
        cmask = falses(countlines(src * ".bim"))
        cmask[colinds] .= true
    end
    m, n = count(rmask), count(cmask)
    bfsrc = SnpArray(src * ".bed")
    # write filtered bed file
    open(des * ".bed", "w+") do io
        write(io, 0x1b6c)
        write(io, 0x01)
        write(io, Matrix{UInt8}(undef, (m + 3) >> 2, n))
    end
    bfdes = SnpArray(des * ".bed", m, "r+")
    bfdes .= @view bfsrc[rmask, cmask]
    # write filtered fam file
    open(des * ".fam", "w") do io
        for (i, line) in enumerate(eachline(src * ".fam"))
            rmask[i] && println(io, line)
        end
    end
    # write filtered bim file
    open(des * ".bim", "w") do io
        for (j, line) in enumerate(eachline(src * ".bim"))
            cmask[j] && println(io, line)
        end
    end
    # output SnpArray
    bfdes
end
