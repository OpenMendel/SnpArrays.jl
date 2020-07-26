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

function _snparray_ax_additive!(out, s::Matrix{UInt8}, v)
    packedstride = size(s, 1)
    @avx for j ∈ eachindex(v)
        for i ∈ eachindex(out)
            k = 2 * ((i-1) & 3)
            block = s[(j-1) * packedstride + ((i-1) >> 2) + 1]
            Aij = (block >> k) & 3
            out[i] += ((Aij >= 2) + (Aij >= 3)) * v[j]
        end
    end
    out
end

function _snparray_ax_dominant!(out, s::Matrix{UInt8}, v)
    packedstride = size(s, 1)
    @avx for j ∈ eachindex(v)
        for i ∈ eachindex(out)
            k = 2 * ((i-1) & 3)
            block = s[(j-1) * packedstride + ((i-1) >> 2) + 1]
            Aij = (block >> k) & 3
            out[i] += (Aij >= 2) * v[j]
        end
    end
    out
end

function _snparray_ax_recessive!(out, s::Matrix{UInt8}, v)
    packedstride = size(s, 1)
    @avx for j ∈ eachindex(v)
        for i ∈ eachindex(out)
            k = 2 * ((i-1) & 3)
            block = s[(j-1) * packedstride + ((i-1) >> 2) + 1]
            Aij = (block >> k) & 3
            out[i] += (Aij >= 3) * v[j]
        end
    end
    out
end

function _snparray_atx_additive!(out, s::Matrix{UInt8}, v)
    packedstride = size(s, 1)
    @avx for i ∈ eachindex(out)
        for j ∈ eachindex(v)
            k = 2 * ((j-1) & 3)
            block = s[(i-1) * packedstride + ((j-1) >> 2) + 1]
            Aij = (block >> k) & 3
            out[i] += ((Aij >= 2) + (Aij >= 3)) * v[j]
        end
    end
    out
end

function _snparray_atx_dominant!(out, s::Matrix{UInt8}, v)
    packedstride = size(s, 1)
    @avx for i ∈ eachindex(out)
        for j ∈ eachindex(v)
            k = 2 * ((j-1) & 3)
            block = s[(i-1) * packedstride + ((j-1) >> 2) + 1]
            Aij = (block >> k) & 3
            out[i] += (Aij >= 2) * v[j]
        end
    end
    out
end

function _snparray_atx_recessive!(out, s::Matrix{UInt8}, v)
    packedstride = size(s, 1)
    @avx for i ∈ eachindex(out)
        for j ∈ eachindex(v)
            k = 2 * ((j-1) & 3)
            block = s[(i-1) * packedstride + ((j-1) >> 2) + 1]
            Aij = (block >> k) & 3
            out[i] += (Aij >= 3) * v[j]
        end
    end
    out
end

function mul!(
    out::AbstractVector{T}, 
    sla::SnpLinAlg{T}, 
    v::AbstractVector{T}) where T <: AbstractFloat
    @assert length(out) == size(sla, 1) && length(v) == size(sla, 2)
    if sla.scale
        sla.storagev2 .= sla.σinv .* v
        w = sla.storagev2
    else
        w = v
    end
    fill!(out, zero(eltype(out)))
    s = sla.s
    if sla.model == ADDITIVE_MODEL
        _snparray_ax_additive!(out, s.data, w)
    elseif sla.model == DOMINANT_MODEL
        _snparray_ax_dominant!(out, s.data, w)
    else
        _snparray_ax_recessive!(out, s.data, w)
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
    @assert length(out) == size(st, 1) && length(v) == size(st, 2)
    sla = st.parent
    s = sla.s
    fill!(out, zero(eltype(out)))
    if sla.model == ADDITIVE_MODEL
        _snparray_atx_additive!(out, s.data, v)
    elseif sla.model == DOMINANT_MODEL
        _snparray_atx_dominant!(out, s.data, v)
    else
        _snparray_atx_recessive!(out, s.data, v)
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
