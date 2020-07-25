struct SnpLinAlg{T} <: AbstractMatrix{T}
    s::SnpArray
    model::Union{Val{1}, Val{2}, Val{3}}
    center::Bool
    scale::Bool
    μ::Vector{T}
    σinv::Vector{T}
    storagev1::Vector{T}
    storagev2::Vector{T}
end

function SnpLinAlg{T}(
    s::AbstractSnpArray;
    model = ADDITIVE_MODEL,
    center::Bool = false,
    scale::Bool = false) where T <: AbstractFloat
    if model == ADDITIVE_MODEL
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
    SnpLinAlg{T}(s, model, center, scale, μ, σinv, storagev1, storagev2)
end

Base.size(sla::SnpLinAlg) = size(sla.s)
Base.size(sla::SnpLinAlg, k::Integer) = size(sla.s, k)

eltype(bm::SnpLinAlg) = eltype(bm.μ)

function mul!(
    out::AbstractVector{T}, 
    sla::SnpLinAlg{T}, 
    v::AbstractVector{T}) where T <: AbstractFloat
    if sla.scale
        sla.storagev2 .= sla.σinv .* v
        w = sla.storagev2
    else
        w = v
    end
    fill!(out, zero(eltype(out)))
    s = sla.s
    if sla.model == ADDITIVE_MODEL
        @avx for j ∈ eachindex(v)
            for i ∈ eachindex(out)
                Aij = s[i, j]
                out[i] += (((Aij >= 2) + (Aij >= 3))) * w[j]
            end
        end
    elseif sla.model == DOMINANT_MODEL
        @avx for j ∈ eachindex(v)
            for i ∈ eachindex(out)
                Aij = s[i, j]
                out[i] += (Aij >= 2) * w[j]
            end
        end
    else
        @avx for j ∈ eachindex(v)
            for i ∈ eachindex(out)
                sij = s[i, j]
                out[i] += (sij >= 3) * w[j]
            end
        end
    end   
    if sla.center
        return out .-= dot(sla.μ, w)
    else
        return out
    end
end

function mul!(
    out::AbstractVector{T}, 
    st::Union{Transpose{T, SnpLinAlg{T}}, Adjoint{T, SnpLinAlg{T}}},
    v::AbstractVector{T}) where T <: AbstractFloat
    sla = st.parent
    s = sla.s
    fill!(out, zero(eltype(out)))
    if sla.model == ADDITIVE_MODEL
        @avx for i ∈ eachindex(out)
            outi = zero(eltype(out))
            for j ∈ eachindex(v)
                stji = s[j, i]
                outi += (((stji >= 2) + (stji >= 3))) * v[j]
            end
            out[i] = outi
        end
    elseif sla.model == DOMINANT_MODEL
        @avx for i ∈ eachindex(out)
            outi = zero(eltype(out))
            for j ∈ eachindex(v)
                stji = s[j, i]
                outi += (stji >= 2) * v[j]
            end
            out[i] = outi
        end
    else
        @avx for i ∈ eachindex(out)
            outi = zero(eltype(out))
            for j ∈ eachindex(v)
                stji = s[j, i]
                outi += (stji >= 3) * v[j]
            end
            out[i] = outi
        end
    end   
    if sla.center
        out .-= sum(v) .* sla.μ
    end
    if sla.scale
        return out .*= sla.σinv
    else
        return out
    end
end
