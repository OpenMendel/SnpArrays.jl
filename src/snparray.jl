"""
    SnpArray

Raw .bed file as a shared, memory-mapped Matrix{UInt8}.  The number of rows, `m`
is stored separately because it is not uniquely determined by the size of the `data` field.
"""
struct SnpArray <: AbstractMatrix{UInt8}
    data::Matrix{UInt8}
    columncounts::Matrix{Int}
    rowcounts::Matrix{Int}
    m::Int
end
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
SnpArray(nm::AbstractString, args...; kwargs...) = SnpArray(nm, countlines(string(splitext(nm)[1], ".fam")), args...; kwargs...)

function SnpArray(::UndefInitializer, m::Integer, n::Integer)
    SnpArray(Matrix{UInt8}(undef, ((m + 3) >> 2, n)), zeros(Int, (4, n)), zeros(Int, (4, m)), m)
end

function SnpArray(file::AbstractString, f::SnpArray)
    open(file, "w+") do io
        write(io, 0x1b6c)
        write(io, 0x01)
        write(io, f.data)
    end
    SnpArray(file, f.m, "r+")
end

function SnpArray(file::AbstractString, m::Integer, n::Integer)
    open(file, "w+") do io
        write(io, 0x1b6c)
        write(io, 0x01)
        write(io, fill(0x00, ((m + 3) >> 2, n)))
    end
    SnpArray(file, m, "r+")
end

StatsBase.counts(f::SnpArray; dims=:) = _counts(f, dims)

function _counts(f::SnpArray, dims::Integer)
    if isone(dims)
        cc = f.columncounts
        if all(iszero, cc)
            m, n = size(f)
            @inbounds for j in 1:n
                for i in 1:m
                    cc[f[i, j] + 1, j] += 1            
                end
            end
        end
        return cc
    elseif dims == 2
        rc = f.rowcounts
        if all(iszero, rc)
            m, n = size(f)
            @inbounds for j in 1:n
                for i in 1:m
                    rc[f[i, j] + 1, i] += 1            
                end
            end
        end
        return rc
    else
        throw(ArgumentError("counts(f::SnpArray, dims=k) only defined for k = 1 or 2"))
    end
end

_counts(f::SnpArray, ::Colon) = sum(_counts(f, 1), dims=2)

function Base.getindex(f::SnpArray, i::Int)  # Linear indexing
    d, r = divrem(i - 1, f.m)
    f[r + 1, d + 1]
end

@inline function Base.getindex(f::SnpArray, i::Integer, j::Integer)
    @boundscheck checkbounds(f, i, j)
    ip3 = i + 3
    (f.data[ip3 >> 2, j] >> ((ip3 & 0x03) << 1)) & 0x03
end

function Base.setindex!(f::SnpArray, x::UInt8, i::Int)  # Linear indexing
    d, r = divrem(i - 1, f.m)
    Base.setindex!(f, x, r + 1, d + 1)
end

@inline function Base.setindex!(f::SnpArray, x::UInt8, i::Integer, j::Integer)
    @boundscheck checkbounds(f, i, j)
    ip3 = i + 3
    shft = (ip3 & 0x03) << 1
    mask = ~(0x03 << shft)
    f.data[ip3 >> 2, j] = (f.data[ip3 >> 2, j] & mask) | (x << shft)
    x
end

Base.eltype(f::SnpArray) = UInt8

Base.length(f::SnpArray) = f.m * size(f.data, 2)

Statistics.mean(f::SnpArray; dims=:) = _mean(f, dims)

function _mean(f::SnpArray,  dims::Integer)
    m, n = size(f)
    if isone(dims)
        cc = _counts(f, 1)   # need to use extractor to force evaluation if needed
        means = Matrix{Float64}(undef, (1, n))
        @inbounds for j in 1:n
            means[j] = (cc[3, j] + 2*cc[4, j]) / (cc[1, j] + cc[3, j] + cc[4, j])
        end
        return means
    elseif dims == 2
        rc = _counts(f, 2)
        means = Matrix{Float64}(undef, (m, 1))
        @inbounds for i in 1:m
            means[i] = (rc[3, i] + 2*rc[4, i]) / (rc[1, i] + rc[3, i] + rc[4, i])
        end
        return means
    else
        throw(ArgumentError("mean(f::SnpArray, dims=k) only defined for k = 1 or 2"))
    end
end

function _mean(f::SnpArray, ::Colon)
    rc = _counts(f, 2)
    (sum(view(rc, 3, :)) + 2*sum(view(rc, 4, :))) / sum(view(rc, [1, 3, 4], :))
end

Base.size(f::SnpArray) = f.m, size(f.data, 2)

Base.size(f::SnpArray, k::Integer) = 
k == 1 ? f.m : k == 2 ? size(f.data, 2) : k > 2 ? 1 : error("Dimension out of range")

Statistics.var(f::SnpArray; corrected::Bool=true, mean=nothing, dims=:) = _var(f, corrected, mean, dims)

function _var(f::SnpArray, corrected::Bool, mean, dims::Integer)
    m, n = size(f)
    means = something(mean, Statistics.mean(f, dims=dims))
    if isone(dims)
        cc = _counts(f, 1)
        vars = Matrix{Float64}(undef, (1, n))
        for j in 1:n
            mnj = means[j]
            vars[j] = (abs2(mnj)*cc[1,j] + abs2(1.0 - mnj)*cc[3,j] + abs2(2.0 - mnj)*cc[4,j]) /
            (cc[1,j] + cc[3,j] + cc[4,j] - (corrected ? 1 : 0))
        end
        return vars
    elseif dims == 2
        rc = _counts(f, 2)
        vars = Matrix{Float64}(undef, (m, 1))
        for i in 1:m
            mni = means[i]
            vars[i] = (abs2(mni)*rc[1,i] + abs2(1.0 - mni)*rc[3,i] + abs2(2.0 - mni)*rc[4,i]) /
            (rc[1,i] + rc[3,i] + rc[4,i] - (corrected ? 1 : 0))
        end
        return vars
    end
    throw(ArgumentError("var(f::SnpArray, dims=k) only defined for k = 1 or 2"))
end

function maf!(out::AbstractVector{T}, f::SnpArray) where T <: AbstractFloat
    cc = _counts(f, 1)
    @inbounds for j in 1:size(f, 2)
        out[j] = (cc[3, j] + 2cc[4, j]) / 2(cc[1, j] + cc[3, j] + cc[4, j])
        (out[j] > 0.5) && (out[j] = 1 - out[j])
    end
    out
end
maf(f::SnpArray) = maf!(Vector{Float64}(undef, size(f, 2)), f)

function minorallele!(out::AbstractVector{Bool}, f::SnpArray)
    cc = _counts(f, 1)
    @inbounds for j in 1:size(f, 2)
        out[j] = cc[1, j] > cc[4, j]
    end
    out
end
minorallele(f::SnpArray) = minorallele!(Vector{Bool}(undef, size(f, 2)), f)

"""    
    outer(f::SnpArray, colinds)
    outer(f::SnpArray)

Return the "outer product", `f * f'` using the `Float32[0, NaN, 1, 2]` encoding of `f`    

The `colinds` argument, when given, causes the operation to be performed on that subset
of the columns.
"""
function outer(f::SnpArray, colinds::AbstractVector{<:Integer})
    m = size(f, 1)
    outer!(Symmetric(zeros(Float32, (m, m))), f, colinds)
end    
outer(f::SnpArray) = outer(f, 1:size(f, 2))

function _copyto_additive!(v::AbstractVector{T}, f::SnpArray, j::Integer) where T <: AbstractFloat
    @inbounds for i in 1:f.m
        fij = f[i, j]
        v[i] = iszero(fij) ? zero(T) : isone(fij) ? T(NaN) : fij - 1
    end    
    v
end

function _copyto_dominant!(v::AbstractVector{T}, f::SnpArray, j::Integer) where T <: AbstractFloat
    @inbounds for i in 1:f.m
        fij = f[i, j]
        v[i] = iszero(fij) ? zero(T) : isone(fij) ? T(NaN) : 1
    end
    v
end

function _copyto_recessive!(v::AbstractVector{T}, f::SnpArray, j::Integer) where T <: AbstractFloat
    @inbounds for i in 1:f.m
        fij = f[i, j]
        v[i] = (iszero(fij) || fij == 2) ? zero(T) : isone(fij) ? T(NaN) : 1
    end    
    v
end

function Base.copyto!(
    v::AbstractVector{T}, 
    f::SnpArray, 
    j::Integer; 
    model::Symbol = :additive,
    center::Bool = false,
    scale::Bool = false,
    impute::Bool = false
    ) where T <: AbstractFloat
    if model == :additive
        _copyto_additive!(v, f, j)
    elseif model == :dominant
        _copyto_dominant!(v, f, j)
    elseif model == :recessive
        _copyto_recessive!(v, f, j)
    else
        throw(ArgumentError("model has to be :additive, :dominant, or :recessive; got $model"))
    end
    if center || scale || impute
        cc = _counts(f, 1)
        μ = model == :additive ? (cc[3, j] + 2cc[4, j]) / (cc[1, j] + cc[3, j] + cc[4, j]) : 
            model == :dominant ? (cc[3, j] +  cc[4, j]) / (cc[1, j] + cc[3, j] + cc[4, j]) :
            cc[4, j] / (cc[1, j] + cc[3, j] + cc[4, j])
        σ = model == :additive ? sqrt(μ * (1 - μ / 2)) : sqrt(μ * (1 - μ))
        doscale = scale && (σ > 0)
        @inbounds for i in 1:f.m
            impute && isnan(v[i]) && (v[i] = T(μ))
            center && (v[i] -= μ)
            doscale && (v[i] /= σ)
        end
    end
    v
end

function Base.copyto!(
    v::AbstractMatrix{T}, 
    f::SnpArray, 
    colinds::AbstractVector{<:Integer};
    kwargs...
    ) where T <: AbstractFloat
    for (vj, j) in enumerate(colinds)
        Base.copyto!(view(v, :, vj), f, j; kwargs...)
    end
    v
end

function Base.copyto!(
    v::AbstractMatrix{T}, 
    f::SnpArray, 
    colmask::AbstractVector{Bool};
    kwargs...
    ) where T <: AbstractFloat
    length(colmask) == size(f, 2) || throw(ArgumentError("`length(colmask)` does not match `size(f, 2)`"))
    vj = 1
    for j in 1:length(colmask)
        if colmask[j] 
            Base.copyto!(view(v, :, vj), f, j; kwargs...)
            vj += 1
        end
    end
    v
end

function Base.convert(t::Type{Vector{T}}, f::SnpArray, j::Integer; kwargs...) where T <: AbstractFloat
    Base.copyto!(Vector{T}(undef, f.m), f, j; kwargs...)
end
function Base.convert(t::Type{Matrix{T}}, f::SnpArray, colinds::AbstractVector{<:Integer}; kwargs...) where T <: AbstractFloat
    Base.copyto!(Matrix{T}(undef, f.m, length(colinds)), f, colinds; kwargs...)
end
function Base.convert(t::Type{Matrix{T}}, f::SnpArray, colmask::AbstractVector{Bool}; kwargs...) where T <: AbstractFloat
    Base.copyto!(Matrix{T}(undef, f.m, count(colmask)), f, colmask; kwargs...)
end
Base.convert(t::Type{Matrix{T}}, f::SnpArray; kwargs...) where T <: AbstractFloat = Base.convert(t, f, 1:size(f, 2); kwargs...)

"""
    outer!(sy::Symmetric, f::SnpArray, colinds)

update `sy` with the sum of the outer products of the columns in `colind` from `f`    
"""
function outer!(sy::Symmetric{T}, f::SnpArray, colinds::AbstractVector{<:Integer}) where T
    tempv = Vector{T}(undef, f.m)
    for j in colinds
        LinearAlgebra.BLAS.syr!(sy.uplo, one(T), copyto!(tempv, f, j), sy.data)
    end    
    sy
end

"""
    missingpos(f::SnpArray)

Return a `SparseMatrixCSC{Bool,Int32}` of the same size as `f` indicating the positions with missing data
"""
function missingpos(f::SnpArray)
    m, n = size(f)
    colptr = sizehint!(Int32[1], n + 1)
    rowval = Int32[]
    @inbounds for j in 1:n
        msngpos = Int32[]
        for i in 1:m
            isone(f[i, j]) && push!(msngpos, i)
        end
        append!(rowval, msngpos)
        push!(colptr, colptr[end] + length(msngpos))
    end
    SparseMatrixCSC(m, n, colptr, rowval, fill(true, length(rowval)))
end

function _missingrate!(out::AbstractVector{<:AbstractFloat}, f::SnpArray,  dims::Integer)
    m, n = size(f)
    if isone(dims)
        cc = _counts(f, 1)   # need to use extractor to force evaluation if needed
        @inbounds for j in 1:n
            out[j] = cc[2, j] / m
        end
    elseif dims == 2
        rc = _counts(f, 2)
        @inbounds for i in 1:m
            out[i] = rc[2, i] / n
        end
    else
        throw(ArgumentError("_missingrate(out, f::SnpArray, dims=k) only defined for k = 1 or 2"))
    end
    out
end
function missingrate(f::SnpArray,  dims::Integer)
    if isone(dims)
        return _missingrate!(Vector{Float64}(undef, size(f, 2)), f, 1)
    elseif dims == 2 
        return _missingrate!(Vector{Float64}(undef, f.m), f, 2)
    else
        throw(ArgumentError("missingrate(f::SnpArray, dims=k) only defined for k = 1 or 2"))
    end
end


"""
    grm(A; method=:GRM, maf_threshold=0.01)

Compute empirical kinship matrix from a SnpArray. Missing genotypes are imputed
on the fly by mean.

# Input  
- `f`: a SnpArray

# Optional Arguments
- `method`: `:GRM` (default), `:MoM`, or `Robust`
- `maf_threshold`: columns with MAF `<maf_threshold` are excluded; default 0.01
- `cinds`: indices or mask of columns to be used for calculating GRM
- `t`: Float type for calculating GRM
"""
function grm(
    f::SnpArray;
    method::Symbol = :GRM,
    maf_threshold::Real = 0.01,
    cinds::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    t::Type{T} = Float64
    ) where T <: AbstractFloat
    mf = maf(f)
    colinds = something(cinds, mf .≥ maf_threshold)
    n = eltype(colinds) == Bool ? count(colinds) : length(colinds)
    G = Mmap.mmap(Matrix{t}, f.m, n)
    if method == :GRM
        Base.copyto!(G, f, colinds, model=:additive, impute=true, center=true, scale=true)
        Φ = G * transpose(G)
        Φ ./= 2n
    elseif method == :MoM
        Base.copyto!(G, f, colinds, model=:additive, impute=true)
        G .-= 1
        Φ = G * transpose(G)
        c = sum(x -> abs2(x) + abs2(1 - x), mf)
        shft, scal = n / 2 - c, 1 / (n - c)
        @inbounds @simd for i in eachindex(Φ)
            Φ[i] = (Φ[i] / 2 + shft) * scal
        end
    elseif method == :Robust
        Base.copyto!(G, f, colinds, model=:additive, center=true, impute=true)
        scal = sum(x -> 4x * (1 - x), mf)
        Φ = G * transpose(G)
        Φ ./= scal
    else
        throw(ArgumentError("method should be :GRM, :MoM, or :Robust; got $method"))
    end
    Φ
end # function grm

"""
    SnpArrays.filter(A[, min_success_rate_per_row, min_success_rate_per_col, maxiters])

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
    f::SnpArray, 
    min_success_rate_per_row::Real = 0.98,
    min_success_rate_per_col::Real = 0.98,
    maxiters::Integer = 5)
    m, n = size(f)
    rc, cc = zeros(Int, m), zeros(Int, n)
    rmask, cmask = trues(m), trues(n)
    rmiss, cmiss = 1 - min_success_rate_per_row, 1 - min_success_rate_per_col
    for iter in 1:maxiters
        fill!(rc, 0)
        fill!(cc, 0)
        @inbounds for j in 1:n
            cmask[j] || continue
            for i in 1:m
                rmask[i] && f[i, j] == 0x01 && (rc[i] += 1; cc[j] += 1)
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
