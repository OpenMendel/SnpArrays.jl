# taken from Gaius.jl
macro _spawn(ex)
    if Threads.nthreads() > 1
        esc(Expr(:macrocall, Expr(:(.), :Threads, QuoteNode(Symbol("@spawn"))), __source__, ex))
    else
        esc(ex)
    end
end
macro _sync(ex)
    if Threads.nthreads() > 1
        esc(Expr(:macrocall, Symbol("@sync"), __source__, ex))
    else
        esc(ex)
    end    
end

struct SnpLinAlg{T} <: AbstractMatrix{T}
    s::SnpArray
    model::Union{Val{1}, Val{2}, Val{3}}
    center::Bool
    scale::Bool
    impute::Bool
    μ::Vector{T}
    σinv::Vector{T}
    storagev1::Vector{T}
    storagev2::Vector{T}
end

"""
    SnpLinAlg{T}(s; model=ADDITIVE_MODEL, center=false, scale=false, impute=true)

Pad a `SnpArray` with some parameters for linear algebera operation.

# Arguments
- s: a `SnpArray`.
- model: one of `ADDITIVE_MODEL`(default), `DOMINANT_MODEL`, `RECESSIVE_MODEL`.
- center: whether to center (default: false).
- scale: whether to scale to standard deviation 1 (default: false).
- impute: whether to impute missing value with column mean (default: true).
"""
function SnpLinAlg{T}(
    s::AbstractSnpArray;
    model = ADDITIVE_MODEL,
    center::Bool = false,
    scale::Bool = false,
    impute::Bool = true) where T <: AbstractFloat
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
    SnpLinAlg{T}(s, model, center, scale, impute, μ, σinv, storagev1, storagev2)
end

Base.size(sla::SnpLinAlg) = size(sla.s)
Base.size(sla::SnpLinAlg, k::Integer) = size(sla.s, k)

eltype(bm::SnpLinAlg) = eltype(bm.μ)


function _snparray_ax_additive_rem!(out, s, v)
    maxp = length(out)
    @avx for j in eachindex(v)
        block = s[1, j]
        for p in 1:maxp
            Aij = (block >> (2*(p-1))) & 3
            out[p] += ((Aij >= 2) + (Aij == 3)) * v[j]
        end
    end
end

function _snparray_ax_additive!(out, s, v, rows, cols, μ)
    #fill!(out, zero(Float64))
    k = rows >> 2
    rem = rows & 3
    
    @avx for j ∈ 1:cols
        for l in 1:k
            block = s[l, j]
            
            for p in 1:4
                Aij = (block >> (2 *(p-1))) & 3
                out[4*(l-1)+p] += ((Aij >= 2) + (Aij == 3)) * v[j]
            end

        end
    end
    if rem != 0
        _snparray_ax_additive_rem!(@view(out[4k+1:end]), @view(s[k+1:k+1, :]), v)
    end
    nothing
end

function _snparray_ax_dominant_rem!(out, s, v)
    maxp = length(out)
    @avx for j in eachindex(v)
        block = s[1, j]
        for p in 1:maxp
            Aij = (block >> (2*(p-1))) & 3
            out[p] += (Aij >= 2) * v[j]
        end
    end
end

function _snparray_ax_dominant!(out, s, v, rows, cols, μ)
    #fill!(out, zero(Float64))
    k = rows >> 2
    rem = rows & 3
    
    @avx for j ∈ 1:cols
        for l in 1:k
            block = s[l, j]
            
            for p in 1:4
                Aij = (block >> (2 *(p-1))) & 3
                out[4*(l-1)+p] += (Aij >= 2) * v[j]
            end

        end
    end
    if rem != 0
        _snparray_ax_dominant_rem!(@view(out[4k+1:end]), @view(s[k+1:k+1, :]), v)
    end
    nothing
end

function _snparray_ax_recessive_rem!(out, s, v)
    maxp = length(out)
    @avx for j in eachindex(v)
        block = s[1, j]
        for p in 1:maxp
            Aij = (block >> (2*(p-1))) & 3
            out[p] += (Aij == 3) * v[j]
        end
    end
end

function _snparray_ax_recessive!(out, s, v, rows, cols, μ)
    #fill!(out, zero(Float64))
    k = rows >> 2
    rem = rows & 3
    
    @avx for j ∈ 1:cols
        for l in 1:k
            block = s[l, j]
            
            for p in 1:4
                Aij = (block >> (2 *(p-1))) & 3
                out[4*(l-1)+p] += (Aij == 3) * v[j]
            end

        end
    end
    if rem != 0
        _snparray_ax_recessive_rem!(@view(out[4k+1:end]), @view(s[k+1:k+1, :]), v)
    end
    nothing
end

function _snparray_ax_additive_meanimpute_rem!(out, s, v, μ)
    maxp = length(out)
    @avx for j in eachindex(v)
        block = s[1, j]
        for p in 1:maxp
            Aij = (block >> (2*(p-1))) & 3
            out[p] += ((Aij >= 2) + (Aij == 3) + (Aij == 1) * μ[j]) * v[j]
        end
    end
end

function _snparray_ax_additive_meanimpute!(out, s, v, rows, cols, μ)
    #fill!(out, zero(Float64))
    k = rows >> 2
    rem = rows & 3
    
    @avx for j ∈ 1:cols
        for l in 1:k
            block = s[l, j]
            
            for p in 1:4
                Aij = (block >> (2 *(p-1))) & 3
                out[4*(l-1)+p] += ((Aij >= 2) + (Aij == 3) + (Aij == 1) * μ[j]) * v[j]
            end

        end
    end
    if rem != 0
        _snparray_ax_additive_meanimpute_rem!(@view(out[4k+1:end]), @view(s[k+1:k+1, :]), v, μ)
    end
    nothing
end

function _snparray_ax_dominant_meanimpute_rem!(out, s, v, μ)
    maxp = length(out)
    @avx for j in eachindex(v)
        block = s[1, j]
        for p in 1:maxp
            Aij = (block >> (2*(p-1))) & 3
            out[p] += ((Aij >= 2) + (Aij == 1) * μ[j]) * v[j]
        end
    end
end

function _snparray_ax_dominant_meanimpute!(out, s, v, rows, cols, μ)
    #fill!(out, zero(Float64))
    k = rows >> 2
    rem = rows & 3
    
    @avx for j ∈ 1:cols
        for l in 1:k
            block = s[l, j]
            
            for p in 1:4
                Aij = (block >> (2 *(p-1))) & 3
                out[4*(l-1)+p] += ((Aij >= 2) + (Aij == 1) * μ[j]) * v[j]
            end

        end
    end
    if rem != 0
        _snparray_ax_dominant_meanimpute_rem!(@view(out[4k+1:end]), @view(s[k+1:k+1, :]), v, μ)
    end
    nothing
end

function _snparray_ax_recessive_meanimpute_rem!(out, s, v, μ)
    maxp = length(out)
    @avx for j in eachindex(v)
        block = s[1, j]
        for p in 1:maxp
            Aij = (block >> (2*(p-1))) & 3
            out[p] += ((Aij == 3) + (Aij == 1) * μ[j]) * v[j]
        end
    end
end

function _snparray_ax_recessive_meanimpute!(out, s, v, rows, cols, μ)
    #fill!(out, zero(Float64))
    k = rows >> 2
    rem = rows & 3
    
    @avx for j ∈ 1:cols
        for l in 1:k
            block = s[l, j]
            
            for p in 1:4
                Aij = (block >> (2 *(p-1))) & 3
                out[4*(l-1)+p] += ((Aij == 3) + (Aij == 1) * μ[j]) * v[j]
            end

        end
    end
    if rem != 0
        _snparray_ax_recessive_meanimpute_rem!(@view(out[4k+1:end]), @view(s[k+1:k+1, :]), v, μ)
    end
    nothing
end

wait(::Nothing) = nothing
function _snparray_ax_tile!(c, A, b, model, μ, impute)
    vstep = 1024
    hstep = 1024
    vstep_log2 = 10
    hstep_log2 = 10

    if !impute
        if model == ADDITIVE_MODEL
            _ftn! = _snparray_ax_additive!
        elseif model == DOMINANT_MODEL
            _ftn! = _snparray_ax_dominant!
        else
            _ftn! = _snparray_ax_recessive!
        end
    else
        if model == ADDITIVE_MODEL
            _ftn! = _snparray_ax_additive_meanimpute!
        elseif model == DOMINANT_MODEL
            _ftn! = _snparray_ax_dominant_meanimpute!
        else
            _ftn! = _snparray_ax_recessive_meanimpute!
        end
    end
    
    fill!(c, zero(UInt8))

    M = length(c) >> 2
    N = size(A, 2)
    Miter = M >>> vstep_log2 # fast div(M, 512)
    Mrem = M & (vstep-1) # fast rem(M, 512)
    Niter = N >>> hstep_log2
    Nrem = N & (hstep-1)
    taskarray = Array{Any}(undef, Miter+1)
    fill!(taskarray, nothing)
    @_sync begin
        GC.@preserve c A b for n in 0:Niter-1
            for m in 0:Miter-1
                wait(taskarray[m+1])
                taskarray[m+1] = @_spawn _ftn!(
                    gesp(stridedpointer(c), ((4 * vstep)*m,)),
                    gesp(stridedpointer(A), (vstep*m, hstep*n)),
                    gesp(stridedpointer(b), (hstep*n,)),
                    vstep << 2,
                    hstep, μ
                )
            end
            if Mrem != 0
                wait(taskarray[Miter+1])
                taskarray[Miter+1] = @_spawn _ftn!(
                    @view(c[(4*vstep)*Miter+1:end]), 
                    @view(A[vstep*(Miter)+1:end, hstep*n+1:hstep*(n+1)]),
                    @view(b[hstep*n+1:hstep*(n+1)]),
                    length(c) - 4*vstep*Miter,
                    hstep, μ
                )
            end
        end
        if Nrem != 0
            for m in 0:Miter-1
                wait(taskarray[m+1])
                taskarray[m+1] = @_spawn _ftn!(
                    @view(c[4*vstep*m+1:4*vstep*(m+1)]),
                    @view(A[vstep*m+1:vstep*(m+1), hstep*Niter+1:end]),
                    @view(b[hstep*Niter+1:end]),
                    vstep << 2, 
                    Nrem, μ
                )
            end
            if Mrem != 0
                wait(taskarray[Miter+1])
                taskarray[Miter+1] = @_spawn _ftn!(
                    @view(c[(4*vstep*Miter+1:end)]),
                    @view(A[vstep*Miter+1:end, hstep*Niter+1:end]),
                    @view(b[hstep*Niter+1:end]),
                    length(c) - 4*vstep*Miter,
                    Nrem, μ
                )
            end
            # s = @_spawn _ftn!(c, @view(A[:, (hstep*Niter+1):end]), @view(b[hstep*Niter+1:end]),
            #     length(c), Nrem, μ
            # )
        end
    end
end

function _snparray_atx_additive_rem!(out, s, v, μ)
    maxp = length(v)
    @avx for i in eachindex(out)
        block = s[1, i]
        for p in 1:maxp
            Aij = (block >> (2*(p-1))) & 3
            out[i] += ((Aij >= 2) + (Aij == 3)) * v[p]
        end
    end
end

function _snparray_atx_additive!(out, s, v, rows, cols, μ)
    #fill!(out, zero(Float64))
    k = rows >> 2
    rem = rows & 3
    
    @avx for i ∈ 1:cols
        for l in 1:k
            block = s[l, i]
            
            for p in 1:4
                Aij = (block >> (2 *(p-1))) & 3
                out[i] += ((Aij >= 2) + (Aij == 3)) * v[4*(l-1)+p]
            end

        end
    end
    if rem != 0
        _snparray_atx_additive_rem!(out, @view(s[k+1:k+1, :]), @view(v[4k+1:end]), μ)
    end
    nothing
end

function _snparray_atx_dominant_rem!(out, s, v, μ)
    maxp = length(v)
    @avx for i in eachindex(out)
        block = s[1, i]
        for p in 1:maxp
            Aij = (block >> (2*(p-1))) & 3
            out[i] += (Aij >= 2) * v[p]
        end
    end
end

function _snparray_atx_dominant!(out, s, v, rows, cols, μ)
    #fill!(out, zero(Float64))
    k = rows >> 2
    rem = rows & 3
    
    @avx for i ∈ 1:cols
        for l in 1:k
            block = s[l, i]
            
            for p in 1:4
                Aij = (block >> (2 *(p-1))) & 3
                out[i] += (Aij >= 2) * v[4*(l-1)+p]
            end

        end
    end
    if rem != 0
        _snparray_atx_dominant_rem!(out, @view(s[k+1:k+1, :]), @view(v[4k+1:end]), μ)
    end
    nothing
end

function _snparray_atx_recessive_rem!(out, s, v, μ)
    maxp = length(v)
    @avx for i in eachindex(out)
        block = s[1, i]
        for p in 1:maxp
            Aij = (block >> (2*(p-1))) & 3
            out[i] += (Aij == 3) * v[p]
        end
    end
end

function _snparray_atx_recessive!(out, s, v, rows, cols, μ)
    #fill!(out, zero(Float64))
    k = rows >> 2
    rem = rows & 3
    
    @avx for i ∈ 1:cols
        for l in 1:k
            block = s[l, i]
            
            for p in 1:4
                Aij = (block >> (2 *(p-1))) & 3
                out[i] += (Aij == 3) * v[4*(l-1)+p]
            end

        end

    end
    if rem != 0
        _snparray_atx_recessive_rem!(out, @view(s[k+1:k+1, :]), @view(v[4k+1:end]), μ)
    end
    nothing
end

function _snparray_atx_additive_meanimpute_rem!(out, s, v, μ)
    maxp = length(v)
    @avx for i in eachindex(out)
        block = s[1, i]
        for p in 1:maxp
            Aij = (block >> (2*(p-1))) & 3
            out[i] += ((Aij >= 2) + (Aij == 3) + (Aij == 1) * μ[i]) * v[p]
        end
    end
end

function _snparray_atx_additive_meanimpute!(out, s, v, rows, cols, μ)
    #fill!(out, zero(Float64))
    k = rows >> 2
    rem = rows & 3
    
    @avx for i ∈ 1:cols
        for l in 1:k
            block = s[l, i]
            
            for p in 1:4
                Aij = (block >> (2 *(p-1))) & 3
                out[i] += ((Aij >= 2) + (Aij == 3) + (Aij == 1) * μ[i]) * v[4*(l-1)+p]
            end

        end
    end
    if rem != 0
        _snparray_atx_additive_meanimpute_rem!(out, @view(s[k+1:k+1, :]), @view(v[4k+1:end]), μ)
    end
    nothing
end

function _snparray_atx_dominant_meanimpute_rem!(out, s, v, μ)
    maxp = length(v)
    @avx for i in eachindex(out)
        block = s[1, i]
        for p in 1:maxp
            Aij = (block >> (2*(p-1))) & 3
            out[i] += ((Aij >= 2) + (Aij == 1) * μ[i]) * v[p]
        end
    end
end

function _snparray_atx_dominant_meanimpute!(out, s, v, rows, cols, μ)
    #fill!(out, zero(Float64))
    k = rows >> 2
    rem = rows & 3
    
    @avx for i ∈ 1:cols
        for l in 1:k
            block = s[l, i]
            
            for p in 1:4
                Aij = (block >> (2 *(p-1))) & 3
                out[i] += ((Aij >= 2) + (Aij == 1) * μ[i]) * v[4*(l-1)+p]
            end

        end
    end
    if rem != 0
        _snparray_atx_dominant_meanimpute_rem!(out, @view(s[k+1:k+1, :]), @view(v[4k+1:end]), μ)
    end
    nothing
end

function _snparray_atx_recessive_meanimpute_rem!(out, s, v, μ)
    maxp = length(v)
    @avx for i in eachindex(out)
        block = s[1, i]
        for p in 1:maxp
            Aij = (block >> (2*(p-1))) & 3
            out[i] += ((Aij == 3) + (Aij == 1) * μ[i]) * v[p]
        end
    end
end

function _snparray_atx_recessive_meanimpute!(out, s, v, rows, cols, μ)
    #fill!(out, zero(Float64))
    k = rows >> 2
    rem = rows & 3
    
    @avx for i ∈ 1:cols
        for l in 1:k
            block = s[l, i]
            
            for p in 1:4
                Aij = (block >> (2 *(p-1))) & 3
                out[i] += ((Aij == 3) + (Aij == 1) * μ[i]) * v[4*(l-1)+p]
            end

        end
    end
    if rem != 0
        _snparray_atx_recessive_meanimpute_rem!(out, @view(s[k+1:k+1, :]), @view(v[4k+1:end]), μ)
    end
    nothing
end

function _snparray_atx_tile!(c, A, b, model, μ, impute)
    vstep = 2048
    hstep = 2048
    vstep_log2 = 11
    hstep_log2 = 11

    if !impute
        if model == ADDITIVE_MODEL
            _ftn! = _snparray_atx_additive!
        elseif model == DOMINANT_MODEL
            _ftn! = _snparray_atx_dominant!
        else
            _ftn! = _snparray_atx_recessive!
        end
    else
        if model == ADDITIVE_MODEL
            _ftn! = _snparray_atx_additive_meanimpute!
        elseif model == DOMINANT_MODEL
            _ftn! = _snparray_atx_dominant_meanimpute!
        else
            _ftn! = _snparray_atx_recessive_meanimpute!
        end
    end
    
    fill!(c, zero(UInt8))

    M = length(b) >> 2
    N = size(A, 2)
    Miter = M >>> vstep_log2 # fast div(M, 512)
    Mrem = M & (vstep-1) # fast rem(M, 512)
    Niter = N >>> hstep_log2
    Nrem = N & (hstep-1)

    taskarray = Array{Any}(undef, Niter+1)
    fill!(taskarray, nothing)
    @_sync begin
        GC.@preserve c A b for m in 0:Miter-1
            for n in 0:Niter-1
                wait(taskarray[n+1])
                @_spawn _ftn!(
                    gesp(stridedpointer(c), (hstep*n,)),
                    gesp(stridedpointer(A), (vstep*m, hstep*n)),
                    gesp(stridedpointer(b), ((4 * vstep)*m,)),
                    vstep << 2,
                    hstep, μ
                )
            end
            if Nrem != 0
                wait(taskarray[Niter+1])
                @_spawn _ftn!(
                    @view(c[hstep*Niter+1:end]),
                    @view(A[vstep*m+1:vstep*(m+1), hstep*Niter+1:end]),
                    @view(b[4*vstep*m+1:4*vstep*(m+1)]),
                    vstep << 2, Nrem, μ
                    )
                # @_spawn _ftn!(@view(c[hstep*Niter+1:end]), @view(A[:, (hstep*Niter+1):end]), b,
                #     length(b), Nrem, μ
                # ) 
            end

        end
        if Mrem != 0
            for n in 0:Niter - 1
                wait(taskarray[n+1])
                @_spawn _ftn!(
                    @view(c[hstep*n+1:hstep*(n+1)]),
                    @view(A[vstep*Miter+1:end, hstep*n+1:hstep*(n+1)]),
                    @view(b[4*vstep*Miter+1:end]),
                    length(b) - 4*vstep*Miter,
                    hstep, μ
                )
            end
            if Nrem != 0
                wait(taskarray[Niter+1])
                @_spawn _ftn!(
                    @view(c[hstep*Niter+1:end]),
                    @view(A[vstep*Miter+1:end, hstep*Niter+1:end]),
                    @view(b[4*vstep*Miter+1:end]),
                    length(b) - 4*vstep*Miter,
                    Nrem, μ
                )
            end
            # @_spawn _ftn!(@view(c[hstep*n+1:hstep*(n+1)]), 
            #     @view(A[vstep*(Miter)+1:end, hstep*n+1:hstep*(n+1)]),
            #     @view(b[(4*vstep)*Miter+1:end]),
            #     length(b) - 4*vstep*Miter,
            #     hstep, μ
            # )
        end
    end
end

"""
    LinearAlgebra.mul!(out, s::SnpLiinAlg, v)

In-place matrix-vector multiplication.
"""
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

    _snparray_ax_tile!(out, s.data, w, sla.model, sla.μ, sla.impute)

    if sla.center
        return out .-= dot(sla.μ, w)
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
    st::Union{Transpose{T, SnpLinAlg{T}}, Adjoint{T, SnpLinAlg{T}}},
    v::AbstractVector{T}) where T <: AbstractFloat
    @assert length(out) == size(st, 1) && length(v) == size(st, 2)
    sla = st.parent
    s = sla.s
    fill!(out, zero(eltype(out)))

    _snparray_atx_tile!(out, s.data, v, sla.model, sla.μ, sla.impute)
    if sla.center
        out .-= sum(v) .* sla.μ
    end
    if sla.scale
        return out .*= sla.σinv
    else
        return out
    end
end
