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
            μ[:] = mean(s, dims=1, model=ADDITIVE_MODEL)
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
            μ[:] = mean(s, dims=1, model=DOMINANT_MODEL)
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
            μ[:] = mean(s, dims=1, model=RECESSIVE_MODEL)
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

eltype(bm::SnpBitMatrix) = eltype(bm.μ)
issymmetric(bm::SnpBitMatrix) = issymmetric(bm.B2) && issymmetric(bm.B1)

function mul!(
    out::AbstractVector{T}, 
    s::SnpBitMatrix{T}, 
    v::AbstractVector{T}) where T <: AbstractFloat
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

function mul!(
    out::AbstractVector{T}, 
    st::Union{Transpose{T, SnpBitMatrix{T}}, Adjoint{T, SnpBitMatrix{T}}},
    v::AbstractVector{T}) where T <: AbstractFloat
    s = st.parent
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

# function mulfast!(
#     out::AbstractVector{T},
#     A::BitMatrix,
#     v::AbstractVector{T}
#     ) where T <: AbstractFloat
#     fill!(out, 0)
#     @inbounds for j in 1:size(A, 2)
#         vj = v[j]
#         @simd for i in 1:size(A, 1)
#             A[i, j] && (out[i] += vj)
#         end
#     end
#     out
# end
