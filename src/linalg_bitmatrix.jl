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
    Base.depwarn("`SnpBitMatrix{T}` is deprecated and will be removed soon. " *
        "Please use `SnpLinAlg{T}` instead.", nothing)
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

function Base.getindex(s::SnpBitMatrix{T}, i::Int, j::Int) where T
    x = s.model == ADDITIVE_MODEL ?
        T(getindex(s.B1, i, j) + getindex(s.B2, i, j)) : T(getindex(s.B1, i, j))
    s.center && (x -= s.μ[j])
    s.scale && (x *= s.σinv[j])
    return x
end

eltype(bm::SnpBitMatrix) = eltype(bm.μ)
issymmetric(bm::SnpBitMatrix) = issymmetric(bm.B2) && issymmetric(bm.B1)

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
    Base.copyto!(v, s)

Copy SnpBitMatrix `s` to numeric vector or matrix `v`. If `s` is centered/scaled,
`v` will be centered/scaled using precomputed column mean `s.μ` and inverse std
`s.σinv`.
"""
function Base.copyto!(
    v::AbstractVecOrMat{T},
    s::AbstractSnpBitMatrix
    ) where T <: AbstractFloat
    m, n = size(s, 1), size(s, 2)
    @inbounds for j in 1:n
        @simd for i in 1:m
            v[i, j] = s[i, j]
        end
    end
    return v
end

"""
    Base.convert(t, s)

Convert a SnpBitMatrix `s` to a numeric vector or matrix of same shape as `s`.
If `s` is centered/scaled, `v` will be centered/scaled using precomputed column
mean `s.μ` and inverse std `s.σinv`.

# Arguments
- `t::Type{AbstractVecOrMat{T}}`: Vector or matrix type.
"""
Base.convert(::Type{T}, s::AbstractSnpBitMatrix) where T <: Array = T(s)
Array{T,N}(s::AbstractSnpBitMatrix) where {T,N} =
    copyto!(Array{T,N}(undef, size(s)), s)
