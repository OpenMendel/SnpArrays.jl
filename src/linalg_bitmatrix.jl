struct SnpBitMatrix{T} <: AbstractMatrix{T}
    B1::BitMatrix
    B2::BitMatrix
    model::Union{Val{1}, Val{2}, Val{3}}
    center::Bool
    scale::Bool
    μ::Vector{T}
    σinv::Vector{T}
    storagev1::Vector{T}
    storagev2::Vector{T}
end

AbstractSnpBitMatrix = Union{SnpBitMatrix, SubArray{T, 1, SnpBitMatrix{T}}, 
    SubArray{T, 2, SnpBitMatrix{T}}} where T

"""
    SnpBitMatrix{T}(s; model=ADDITIVE_MODEL, center=false, scale=false)

Store an `AbstractSnpArray` `s` as two `BitMatrix`es to simplify linear algebra.

# Arguments
- s: an `AbstractSnpArray`.
- model: one of `ADDITIVE_MODEL`(default), `DOMINANT_MODEL`, `RECESSIVE_MODEL`.
- center: whether to center (default: false).
- scale: whether to scale to standard deviation 1 (default: false).
"""
function SnpBitMatrix{T}(
    s::AbstractSnpArray;
    model = ADDITIVE_MODEL,
    center::Bool = false,
    scale::Bool = false) where T <: AbstractFloat
    if model == ADDITIVE_MODEL
        B1 = s .≥ 0x02
        B2 = s .≥ 0x03
        if center || scale
            μ = Vector{T}(undef, size(s, 2))
            μ[:] = mean(s; dims=1, model=ADDITIVE_MODEL)
        else
            μ = T[]
        end
        if scale
            σinv = Vector{T}(undef, size(s, 2))
            @inbounds @simd for j in 1:size(s, 2)
                σinv[j] = sqrt(μ[j] * (1 - μ[j] / 2))
                σinv[j] = σinv[j] > 0 ? inv(σinv[j]) : one(T)
            end
        else
            σinv = T[]
        end
    elseif model == DOMINANT_MODEL
        B1 = s .≥ 0x02
        B2 = falses(0, 0)
        if center || scale
            μ = Vector{T}(undef, size(s, 2))
            μ[:] = mean(s; dims=1, model=DOMINANT_MODEL)
        else
            μ = T[]
        end
        if scale
            σinv = Vector{T}(undef, size(s, 2))
            @inbounds @simd for j in 1:size(s, 2)
                σinv[j] = sqrt(μ[j] * (1 - μ[j]))
                σinv[j] = σinv[j] > 0 ? inv(σinv[j]) : one(T)
            end
        else
            σinv = T[]
        end
    elseif model == RECESSIVE_MODEL
        B1 = s .== 0x03
        B2 = falses(0, 0)
        if center || scale
            μ = Vector{T}(undef, size(s, 2))
            μ[:] = mean(s; dims=1, model=RECESSIVE_MODEL)
        else
            μ = T[]
        end
        if scale
            σinv = Vector{T}(undef, size(s, 2))
            @inbounds @simd for j in 1:size(s, 2)
                σinv[j] = sqrt(μ[j] * (1 - μ[j]))
                σinv[j] = σinv[j] > 0 ? inv(σinv[j]) : one(T)
            end
        else
            σinv = T[]
        end
    else
        throw(ArgumentError("unrecognized model $model"))
    end
    storagev1 = Vector{T}(undef, size(s, 1))
    storagev2 = Vector{T}(undef, size(s, 2))
    SnpBitMatrix{T}(B1, B2, model, center, scale, μ, σinv, storagev1, storagev2)
end

Base.size(bm::SnpBitMatrix) = size(bm.B1)
Base.size(bm::SnpBitMatrix, k::Integer) = size(bm.B1, k)

function Base.getindex(s::SnpBitMatrix{T}, i::Int) where T 
    s.model == ADDITIVE_MODEL ?
        T(getindex(s.B1, i) + getindex(s.B2, i)) : T(getindex(s.B1, i))
end
function Base.getindex(s::SnpBitMatrix{T}, i::Int, j::Int) where T 
    s.model == ADDITIVE_MODEL ? 
        T(getindex(s.B1, i, j) + getindex(s.B2, i, j)) : T(getindex(s.B1, i, j))
end

eltype(bm::SnpBitMatrix) = eltype(bm.μ)
issymmetric(bm::SnpBitMatrix) = issymmetric(bm.B2) && issymmetric(bm.B1)

function _gemv_avx_sized!(c, A, b, rows, cols)
    @avx for j in 1:cols, i in 1:rows
        c[i] += A[i, j] * b[j]
    end
end

function _gemv_tile!(c, A, b)
    vstep = 512
    hstep = 512
    vstep_log2 = 9
    hstep_log2 = 9
    M, N = size(A)
    Miter = M >>> vstep_log2 # fast div(M, 512)
    Mrem = M & (vstep-1) # fast rem(M, 512)
    Niter = N >>> hstep_log2
    Nrem = N & (hstep-1)
    GC.@preserve c A b for n in 0:Niter-1
        for m in 0:Miter-1
            _gemv_avx_sized!(
                gesp(stridedpointer(c), (vstep*m,)),
                gesp(stridedpointer(A), (vstep*m, hstep*n)),
                gesp(stridedpointer(b), (hstep*n,)),
                vstep,
                hstep
            )
        end
        if Mrem != 0
            _gemv_avx_sized!(
                gesp(stridedpointer(c), (vstep*Miter,)), 
                gesp(stridedpointer(A), (vstep*Miter, hstep*n)),
                gesp(stridedpointer(b), (hstep*n,)),
                Mrem, hstep
            )
        end
    end
    if Nrem != 0
        _gemv_avx_sized!(
            gesp(stridedpointer(c), (0,)),
            gesp(stridedpointer(A), (0, hstep*Niter)),
            gesp(stridedpointer(b), (hstep*Niter,)),
            length(c),
            Nrem
        )
    end
end

function mul!(c::Vector{T}, A::BitMatrix, b::Vector{T}) where T
    fill!(c, zero(eltype(c)))
    _gemv_tile!(c, A, b)
end

function mul!(c::Vector{T}, A::Union{Transpose{Bool,BitArray{2}}, Adjoint{Bool, BitArray{2}}}, b::Vector{T}) where T
    tA = transpose(A)
    fill!(c, zero(eltype(c)))
    @avx for i in 1:size(tA, 2)
        for j in 1:size(tA, 1)
            c[i] += tA[j, i] * b[j]
        end
    end
end

"""
    LinearAlgebra.mul!(out, s::SnpBitMatrix, v)

In-place matrix-vector multiplication.
"""
function mul!(
    out::AbstractVector{T}, 
    s::SnpBitMatrix{T}, 
    v::AbstractVector{T}) where T <: AbstractFloat
    @assert length(out) == size(s, 1) && length(v) == size(s, 2)
    if s.scale
        s.storagev2 .= s.σinv .* v
        w = s.storagev2
    else
        w = v
    end
    if s.model == ADDITIVE_MODEL
        mul!(out, s.B1, w)
        mul!(s.storagev1, s.B2, w)
        out .+= s.storagev1
    else
        mul!(out, s.B1, w)
    end   
    if s.center
        return out .-= dot(s.μ, w)
    else
        return out
    end
end

"""
    LinearAlgebra.mul!(out, s::Union{Transpose{T, SnpBitMatrix{T}}, Adjoint{T, SnpBitMatrix{T}}}, v)

In-place matrix-vector multiplication, with transposed `SnpBitMatrix`.
"""
function mul!(
    out::AbstractVector{T}, 
    st::Union{Transpose{T, SnpBitMatrix{T}}, Adjoint{T, SnpBitMatrix{T}}},
    v::AbstractVector{T}) where T <: AbstractFloat
    s = st.parent
    @assert length(out) == size(st, 1) && length(v) == size(st, 2)
    if s.model == ADDITIVE_MODEL
        mul!(out, transpose(s.B1), v)
        mul!(s.storagev2, transpose(s.B2), v)
        out .+= s.storagev2
    else
        mul!(out, transpose(s.B1), v)
    end   
    if s.center
        out .-= sum(v) .* s.μ
    end
    if s.scale
        return out .*= s.σinv
    else
        return out
    end
end

"""
    Base.copyto!(v, s, center=false, scale=false)

Copy SnpBitMatrix `s` to numeric vector or matrix `v`. 

# Arguments
- `center::Bool=false`: center column by mean.
- `scale::Bool=false`: scale column by theoretical variance.
"""
function Base.copyto!(
    v::AbstractVecOrMat{T}, 
    s::AbstractSnpBitMatrix;
    center::Bool = false,
    scale::Bool = false,
    ) where T <: AbstractFloat
    m, n = size(s, 1), size(s, 2)
    # no center or scale
    if !center && !scale
        @inbounds for j in 1:n
            @simd for i in 1:m
                v[i, j] = convert(T, s[i, j])
            end
        end
        return v
    end
    # center, scale (TODO: how to use precomputed μ and σinv?)
    model = typeof(s) <: SubArray ? s.parent.model : s.model
    @inbounds for j in 1:n
        μj, mj = zero(T), 0
        @simd for i in 1:m
            vij = convert(T, s[i, j])
            v[i, j] = vij
            μj += isnan(vij) ? zero(T) : vij
            mj += isnan(vij) ? 0 : 1
        end
        μj /= mj
        σj = model == ADDITIVE_MODEL ? sqrt(μj * (1 - μj / 2)) : sqrt(μj * (1 - μj))
        @simd for i in 1:m
            center && (v[i, j] -= μj)
            scale && σj > 0 && (v[i, j] /= σj)
        end
    end
    return v
end

"""
    Base.convert(t, s, model=ADDITIVE_MODEL, center=false, scale=false, impute=false)

Convert a SnpBitMatrix `s` to a numeric vector or matrix of same shape as `s`.

# Arguments
- `t::Type{AbstractVecOrMat{T}}`: Vector or matrix type.
- `model::Union{Val{1}, Val{2}, Val{3}}=ADDITIVE_MODEL`: `ADDITIVE_MODEL` (default), `DOMINANT_MODEL`, or `RECESSIVE_MODEL`.  
- `center::Bool=false`: center column by mean.
- `scale::Bool=false`: scale column by theoretical variance.
"""
function Base.convert(
    ::Type{T},
    s::AbstractSnpBitMatrix;
    kwargs...) where T <: Array
    T(s; kwargs...)
end
Array{T,N}(s::AbstractSnpBitMatrix; kwargs...) where {T,N} = 
    copyto!(Array{T,N}(undef, size(s)), s; kwargs...)
