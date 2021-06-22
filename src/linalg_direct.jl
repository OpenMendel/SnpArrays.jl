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

AbstractSnpLinAlg = Union{SnpLinAlg, SubArray{T, 1, SnpLinAlg{T}}, 
    SubArray{T, 2, SnpLinAlg{T}}} where T

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
    μ = dropdims(mean(s; dims=1, model=model), dims=1)
    σinv = Vector{T}(undef, size(s, 2))
    if model == ADDITIVE_MODEL
        @inbounds @simd for j in 1:size(s, 2)
            σinv[j] = sqrt(μ[j] * (1 - μ[j] / 2))
            σinv[j] = σinv[j] > 0 ? inv(σinv[j]) : one(T)
        end
    elseif model == DOMINANT_MODEL || model == RECESSIVE_MODEL
        @inbounds @simd for j in 1:size(s, 2)
            σinv[j] = sqrt(μ[j] * (1 - μ[j]))
            σinv[j] = σinv[j] > 0 ? inv(σinv[j]) : one(T)
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

function Base.getindex(s::SnpLinAlg{T}, i::Int, j::Int) where T
    x = SnpArrays.convert(T, getindex(s.s, i, j), s.model)
    s.impute && isnan(x) && return s.μ[j]
    s.center && (x -= s.μ[j])
    s.scale && (x *= s.σinv[j])
    return x
end

# macros taken from Gaius.jl
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

"""
    LinearAlgebra.mul!(out, sla::SnpLinAlg, v)

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

    _snparray_ax_tile!(out, s.data, w, sla.model, sla.μ, sla.impute, s.m)

    if sla.center
        return out .-= dot(sla.μ, w)
    else
        return out
    end
end

"""
    LinearAlgebra.mul!(out, sla::SnpLinAlg, V)

In-place matrix-matrix multiplication.
"""
function mul!(
    out::AbstractMatrix{T}, 
    sla::SnpLinAlg{T}, 
    V::AbstractMatrix{T}) where T <: AbstractFloat
    @assert size(out, 1) == size(sla, 1) && size(out, 2) == size(V, 2) && size(sla, 2) == size(V, 1)

    sla.storagev2 .= sla.scale ? sla.σinv : one(T)
    fill!(out, zero(eltype(out)))
    s = sla.s

    _snparray_AX_tile!(out, s.data, V, sla.model, sla.μ, sla.impute, s.m, sla.storagev2)

    return out
end

"""
    LinearAlgebra.mul!(out, st::Union{Transpose{T, SnpLinAlg{T}}, Adjoint{T, SnpLinAlg{T}}}, v)

In-place matrix-vector multiplication, with transposed `SnpLinAlg`.
"""
function mul!(
    out::AbstractVector{T}, 
    st::Union{Transpose{T, SnpLinAlg{T}}, Adjoint{T, SnpLinAlg{T}}},
    v::AbstractVector{T}) where T <: AbstractFloat
    @assert length(out) == size(st, 1) && length(v) == size(st, 2)
    sla = st.parent
    s = sla.s
    fill!(out, zero(eltype(out)))

    _snparray_atx_tile!(out, s.data, v, sla.model, sla.μ, sla.impute, s.m)
    if sla.center
        out .-= sum(v) .* sla.μ
    end
    if sla.scale
        return out .*= sla.σinv
    else
        return out
    end
end

"""
    LinearAlgebra.mul!(out, st::Union{Transpose{T, SnpLinAlg{T}}, Adjoint{T, SnpLinAlg{T}}}, V)

In-place matrix-matrix multiplication, with transposed `SnpLinAlg`.
"""
function mul!(
    out::AbstractMatrix{T}, 
    st::Union{Transpose{T, SnpLinAlg{T}}, Adjoint{T, SnpLinAlg{T}}},
    V::AbstractMatrix{T}) where T <: AbstractFloat
    sla = st.parent
    @assert size(out, 1) == size(sla, 2) && size(out, 2) == size(V, 2) && size(sla, 1) == size(V, 1)

    s = sla.s
    fill!(out, zero(eltype(out)))
    sla.storagev2 .= sla.scale ? sla.σinv : one(T)

    _snparray_AtX_tile!(out, s.data, V, sla.model, sla.μ, sla.impute, s.m, sla.storagev2)

    return out
end

wait(::Nothing) = nothing

function _snparray_ax_tile!(c, A, b, model, μ, impute, rows_filled)
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

    M = length(c) >> 2
    N = size(A, 2)
    Miter = M >>> vstep_log2 # fast div(M, 1024)
    Mrem = rows_filled & (vstep << 2 - 1) # fast rem(rows_filled, 4vstep)
    Niter = N >>> hstep_log2
    Nrem = N & (hstep - 1)
    taskarray = Array{Any}(undef, Miter + 1)
    fill!(taskarray, nothing)
    @_sync begin
        GC.@preserve c A b for n in 0:Niter - 1
            for m in 0:Miter - 1
                wait(taskarray[m+1])
                taskarray[m+1] = @_spawn _ftn!(
                    gesp(stridedpointer(c), (4 * vstep * m,)),
                    gesp(stridedpointer(A), (vstep * m, hstep * n)),
                    gesp(stridedpointer(b), (hstep * n,)),
                    vstep << 2, hstep, @view(μ[hstep * n + 1:hstep * (n + 1)])
                )
            end
            if Mrem != 0
                wait(taskarray[Miter+1])
                taskarray[Miter+1] = @_spawn _ftn!(
                    @view(c[4 * vstep * Miter + 1:end]), 
                    @view(A[vstep * Miter + 1:end, hstep * n + 1:hstep * (n + 1)]),
                    @view(b[hstep * n + 1:hstep * (n + 1)]),
                    length(c) - 4 * vstep * Miter, hstep,
                    @view(μ[hstep * n + 1:hstep * (n + 1)])
                )
            end
        end
        if Nrem != 0
            for m in 0:Miter-1
                wait(taskarray[m+1])
                taskarray[m+1] = @_spawn _ftn!(
                    @view(c[4 * vstep * m + 1:4 * vstep * (m + 1)]),
                    @view(A[vstep * m + 1:vstep * (m + 1), hstep * Niter + 1:end]),
                    @view(b[hstep * Niter + 1:end]),
                    vstep << 2, 
                    Nrem, @view(μ[hstep * Niter + 1:end])
                )
            end
            if Mrem != 0
                wait(taskarray[Miter + 1])
                taskarray[Miter + 1] = @_spawn _ftn!(
                    @view(c[4 * vstep * Miter+1:end]),
                    @view(A[vstep * Miter + 1:end, hstep * Niter + 1:end]),
                    @view(b[hstep * Niter + 1:end]),
                    length(c) - 4 * vstep * Miter,
                    Nrem, @view(μ[hstep * Niter + 1:end])
                )
            end
        end
    end
end

function _snparray_AX_tile!(C, A, B, model, μ, impute, rows_filled, σinv)
    vstep = 256
    hstep = 256
    pstep = 256
    vstep_log2 = 8
    hstep_log2 = 8
    pstep_log2 = 8

    if !impute
        if model == ADDITIVE_MODEL
            _ftn! = _snparray_AX_additive!
        elseif model == DOMINANT_MODEL
            _ftn! = _snparray_AX_dominant!
        else
            _ftn! = _snparray_AX_recessive!
        end
    else
        if model == ADDITIVE_MODEL
            _ftn! = _snparray_AX_additive_meanimpute!
        elseif model == DOMINANT_MODEL
            _ftn! = _snparray_AX_dominant_meanimpute!
        else
            _ftn! = _snparray_AX_recessive_meanimpute!
        end
    end

    M = size(C, 1) >> 2
    N = size(A, 2)
    P = size(C, 2)
    Miter = M >>> vstep_log2 # fast div(M, 1024)
    Mrem = rows_filled & (vstep << 2 - 1) # fast rem(rows_filled, 4vstep)
    Niter = N >>> hstep_log2
    Nrem = N & (hstep - 1)
    Piter = P >>> pstep_log2
    Prem = P & (pstep - 1)
    taskarray = Array{Any}(undef, Miter + 1)
    fill!(taskarray, nothing)
    @_sync begin
        GC.@preserve C A B for p in 0:Piter - 1
            for n in 0:Niter - 1
                for m in 0:Miter - 1
                    wait(taskarray[m+1])
                    taskarray[m+1] = @_spawn _ftn!(
                        @view(C[4 * vstep * m + 1:4 * vstep * (m + 1), pstep * p + 1:pstep * (p + 1)]), 
                        @view(A[vstep * m + 1:vstep * (m + 1), hstep * n + 1:hstep * (n + 1)]),
                        @view(B[hstep * n + 1:hstep * (n + 1), pstep * p + 1:pstep * (p + 1)]),
                        vstep << 2, hstep, pstep, 
                        @view(μ[hstep * n + 1:hstep * (n + 1)]),
                        @view(σinv[hstep * n + 1:hstep * (n + 1)])
                    )
                end
                if Mrem != 0
                    wait(taskarray[Miter+1])
                    taskarray[Miter+1] = @_spawn _ftn!(
                        @view(C[4 * vstep * Miter + 1:end, pstep * p + 1:pstep * (p + 1)]), 
                        @view(A[vstep * Miter + 1:end, hstep * n + 1:hstep * (n + 1)]),
                        @view(B[hstep * n + 1:hstep * (n + 1), pstep * p + 1:pstep * (p + 1)]),
                        size(C, 1) - 4 * vstep * Miter, hstep, pstep,
                        @view(μ[hstep * n + 1:hstep * (n + 1)]),
                        @view(σinv[hstep * n + 1:hstep * (n + 1)])
                    )
                end
            end
            if Nrem != 0
                for m in 0:Miter-1
                    wait(taskarray[m+1])
                    taskarray[m+1] = @_spawn _ftn!(
                        @view(C[4 * vstep * m + 1:4 * vstep * (m + 1), pstep * p + 1:pstep * (p + 1)]),
                        @view(A[vstep * m + 1:vstep * (m + 1), hstep * Niter + 1:end]),
                        @view(B[hstep * Niter + 1:end, pstep * p + 1:pstep * (p + 1)]),
                        vstep << 2, Nrem, pstep,
                        @view(μ[hstep * Niter + 1:end]),
                        @view(σinv[hstep * Niter + 1:end])
                    )
                end
                if Mrem != 0
                    wait(taskarray[Miter + 1])
                    taskarray[Miter + 1] = @_spawn _ftn!(
                        @view(C[4 * vstep * Miter+1:end, pstep * p + 1:pstep * (p + 1)]),
                        @view(A[vstep * Miter + 1:end, hstep * Niter + 1:end]),
                        @view(B[hstep * Niter + 1:end, pstep * p + 1:pstep * (p + 1)]),
                        size(C, 1) - 4 * vstep * Miter, Nrem, pstep,
                        @view(μ[hstep * Niter + 1:end]),
                        @view(σinv[hstep * Niter + 1:end])
                    )
                end
            end
        end
        if Prem != 0
            for n in 0:Niter - 1
                for m in 0:Miter - 1
                    wait(taskarray[m+1])
                    taskarray[m+1] = @_spawn _ftn!(
                        @view(C[4 * vstep * m + 1:4 * vstep * (m + 1), pstep * Piter + 1:end]), 
                        @view(A[vstep * m + 1:vstep * (m + 1), hstep * n + 1:hstep * (n + 1)]),
                        @view(B[hstep * n + 1:hstep * (n + 1), pstep * Piter + 1:end]),
                        vstep << 2, hstep, Prem,
                        @view(μ[hstep * n + 1:hstep * (n + 1)]),
                        @view(σinv[hstep * n + 1:hstep * (n + 1)])
                    )
                end
                if Mrem != 0
                    wait(taskarray[Miter+1])
                    taskarray[Miter+1] = @_spawn _ftn!(
                        @view(C[4 * vstep * Miter + 1:end, pstep * Piter + 1:end]), 
                        @view(A[vstep * Miter + 1:end, hstep * n + 1:hstep * (n + 1)]),
                        @view(B[hstep * n + 1:hstep * (n + 1), pstep * Piter + 1:end]),
                        size(C, 1) - 4 * vstep * Miter, hstep, Prem,
                        @view(μ[hstep * n + 1:hstep * (n + 1)]),
                        @view(σinv[hstep * n + 1:hstep * (n + 1)])
                    )
                end
            end
            if Nrem != 0
                for m in 0:Miter-1
                    wait(taskarray[m+1])
                    taskarray[m+1] = @_spawn _ftn!(
                        @view(C[4 * vstep * m + 1:4 * vstep * (m + 1), pstep * Piter + 1:end]),
                        @view(A[vstep * m + 1:vstep * (m + 1), hstep * Niter + 1:end]),
                        @view(B[hstep * Niter + 1:end, pstep * Piter + 1:end]),
                        vstep << 2, Nrem, Prem,
                        @view(μ[hstep * Niter + 1:end]),
                        @view(σinv[hstep * Niter + 1:end])
                    )
                end
                if Mrem != 0
                    wait(taskarray[Miter + 1])
                    taskarray[Miter + 1] = @_spawn _ftn!(
                        @view(C[4 * vstep * Miter+1:end, pstep * Piter + 1:end]),
                        @view(A[vstep * Miter + 1:end, hstep * Niter + 1:end]),
                        @view(B[hstep * Niter + 1:end, pstep * Piter + 1:end]),
                        size(C, 1) - 4 * vstep * Miter, Nrem, Prem,
                        @view(μ[hstep * Niter + 1:end]),
                        @view(σinv[hstep * Niter + 1:end])
                    )
                end
            end
        end
    end
end

function _snparray_atx_tile!(c, A, b, model, μ, impute, rows_filled)
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

    M = length(b) >> 2
    N = size(A, 2)
    Miter = M >>> vstep_log2 # fast div(M, 1024)
    Mrem = rows_filled & (vstep << 2 - 1) # fast rem(rows_filled, 4vstep)
    Niter = N >>> hstep_log2
    Nrem = N & (hstep - 1)
    taskarray = Array{Any}(undef, Niter+1)
    fill!(taskarray, nothing)
    @_sync begin
        GC.@preserve c A b for m in 0:Miter - 1
            for n in 0:Niter - 1
                wait(taskarray[n + 1])
                taskarray[n + 1] = @_spawn _ftn!(
                    gesp(stridedpointer(c), (hstep * n,)),
                    gesp(stridedpointer(A), (vstep * m, hstep * n)),
                    gesp(stridedpointer(b), (4 * vstep * m,)),
                    vstep << 2, hstep, @view(μ[hstep * n + 1:hstep * (n + 1)])
                )
            end
            if Nrem != 0
                wait(taskarray[Niter + 1])
                taskarray[Niter + 1] = @_spawn _ftn!(
                    @view(c[hstep * Niter + 1:end]),
                    @view(A[vstep * m + 1:vstep * (m + 1), hstep * Niter + 1:end]),
                    @view(b[4 * vstep * m + 1:4 * vstep * (m + 1)]),
                    vstep << 2, Nrem, @view(μ[hstep * Niter + 1:end])
                )
            end

        end
        if Mrem != 0
            for n in 0:Niter - 1
                wait(taskarray[n + 1])
                taskarray[n + 1] = @_spawn _ftn!(
                    @view(c[hstep * n + 1:hstep * (n + 1)]),
                    @view(A[vstep * Miter + 1:end, hstep * n + 1:hstep * (n + 1)]),
                    @view(b[4 * vstep * Miter + 1:end]),
                    length(b) - 4 * vstep * Miter,
                    hstep, @view(μ[hstep * n + 1:hstep * (n + 1)])
                )
            end
            if Nrem != 0
                wait(taskarray[Niter + 1])
                taskarray[Niter + 1] = @_spawn _ftn!(
                    @view(c[hstep * Niter + 1:end]),
                    @view(A[vstep * Miter + 1:end, hstep * Niter + 1:end]),
                    @view(b[4 * vstep * Miter + 1:end]),
                    length(b) - 4 * vstep * Miter,
                    Nrem, @view(μ[hstep * Niter + 1:end])
                )
            end
        end
    end
end


function _snparray_AtX_tile!(C, A, B, model, μ, impute, rows_filled, σinv)
    vstep = 2048
    hstep = 2048
    pstep = 2048
    vstep_log2 = 11
    hstep_log2 = 11
    pstep_log2 = 11

    if !impute
        if model == ADDITIVE_MODEL
            _ftn! = _snparray_AtX_additive!
        elseif model == DOMINANT_MODEL
            _ftn! = _snparray_AtX_dominant!
        else
            _ftn! = _snparray_AtX_recessive!
        end
    else
        if model == ADDITIVE_MODEL
            _ftn! = _snparray_AtX_additive_meanimpute!
        elseif model == DOMINANT_MODEL
            _ftn! = _snparray_AtX_dominant_meanimpute!
        else
            _ftn! = _snparray_AtX_recessive_meanimpute!
        end
    end

    M = size(B, 1) >> 2
    N = size(A, 2)
    P = size(C, 2)
    Miter = M >>> vstep_log2 # fast div(M, 1024)
    Mrem = rows_filled & (vstep << 2 - 1) # fast rem(rows_filled, 4vstep)
    Niter = N >>> hstep_log2
    Nrem = N & (hstep - 1)
    Piter = P >>> pstep_log2
    Prem = P & (pstep - 1)
    taskarray = Array{Any}(undef, Niter + 1)
    fill!(taskarray, nothing)
    @_sync begin
        GC.@preserve C A B for p in 0:Piter - 1
            for m in 0:Miter - 1
                for n in 0:Niter - 1
                    wait(taskarray[n + 1])
                    taskarray[n + 1] = @_spawn _ftn!(
                        @view(C[hstep * n + 1:hstep * (n + 1), pstep * p + 1:pstep * (p + 1)]), 
                        @view(A[vstep * m + 1:vstep * (m + 1), hstep * n + 1:hstep * (n + 1)]),
                        @view(B[4 * vstep * m + 1:4 * vstep * (m + 1), pstep * p + 1:pstep * (p + 1)]),
                        vstep << 2, hstep, pstep, 
                        @view(μ[hstep * n + 1:hstep * (n + 1)]),
                        @view(σinv[hstep * n + 1:hstep * (n + 1)])
                    )
                end
                if Nrem != 0
                    wait(taskarray[Niter + 1])
                    taskarray[Niter + 1] = @_spawn _ftn!(
                        @view(C[hstep * Niter + 1:end, pstep * p + 1:pstep * (p + 1)]), 
                        @view(A[vstep * m + 1:vstep * (m + 1), hstep * Niter + 1:end]),
                        @view(B[4 * vstep * m + 1:4 * vstep * (m + 1), pstep * p + 1:pstep * (p + 1)]),
                        vstep << 2, Nrem, pstep,
                        @view(μ[hstep * Niter + 1:end]),
                        @view(σinv[hstep * Niter + 1:end])
                    )
                end
            end
            if Mrem != 0
                for n in 0:Niter - 1
                    wait(taskarray[n + 1])
                    taskarray[n + 1] = @_spawn _ftn!(
                        @view(C[hstep * n + 1:hstep * (n + 1), pstep * p + 1:pstep * (p + 1)]),
                        @view(A[vstep * Miter + 1:end, hstep * n + 1:hstep * (n + 1)]),
                        @view(B[4 * vstep * Miter + 1:end, pstep * p + 1:pstep * (p + 1)]),
                        size(B, 1) - 4 * vstep * Miter, hstep, pstep,
                        @view(μ[hstep * n + 1:hstep * (n + 1)]),
                        @view(σinv[hstep * n + 1:hstep * (n + 1)])
                    )
                end
                if Nrem != 0
                    wait(taskarray[Niter + 1])
                    taskarray[Niter + 1] = @_spawn _ftn!(
                        @view(C[hstep * Niter + 1:end, pstep * p + 1:pstep * (p + 1)]),
                        @view(A[vstep * Miter + 1:end, hstep * Niter + 1:end]),
                        @view(B[4 * vstep * Miter + 1:end, pstep * p + 1:pstep * (p + 1)]),
                        size(B, 1) - 4 * vstep * Miter, Nrem, pstep,
                        @view(μ[hstep * Niter + 1:end]),
                        @view(σinv[hstep * Niter + 1:end])
                    )
                end
            end
        end
        if Prem != 0
            for m in 0:Miter - 1
                for n in 0:Niter - 1
                    wait(taskarray[n + 1])
                    taskarray[n + 1] = @_spawn _ftn!(
                        @view(C[hstep * n + 1:hstep * (n + 1), pstep * Piter + 1:end]), 
                        @view(A[vstep * m + 1:vstep * (m + 1), hstep * n + 1:hstep * (n + 1)]),
                        @view(B[4 * vstep * m + 1:4 * vstep * (m + 1), pstep * Piter + 1:end]),
                        vstep << 2, hstep, Prem,
                        @view(μ[hstep * n + 1:hstep * (n + 1)]),
                        @view(σinv[hstep * n + 1:hstep * (n + 1)])
                    )
                end
                if Nrem != 0
                    wait(taskarray[Niter + 1])
                    taskarray[Niter + 1] = @_spawn _ftn!(
                        @view(C[hstep * Niter + 1:end, pstep * Piter + 1:end]), 
                        @view(A[vstep * m + 1:vstep * (m + 1), hstep * Niter + 1:end]),
                        @view(B[4 * vstep * m + 1:4 * vstep * (m + 1), pstep * Piter + 1:end]),
                        vstep << 2, Nrem, Prem,
                        @view(μ[hstep * Niter + 1:end]),
                        @view(σinv[hstep * Niter + 1:end])
                    )
                end
            end
            if Mrem != 0
                for n in 0:Niter - 1
                    wait(taskarray[n + 1])
                    taskarray[n + 1] = @_spawn _ftn!(
                        @view(C[hstep * n + 1:hstep * (n + 1), pstep * Piter + 1:end]),
                        @view(A[vstep * Miter + 1:end, hstep * n + 1:hstep * (n + 1)]),
                        @view(B[4 * vstep * Miter + 1:end, pstep * Piter + 1:end]),
                        size(B, 1) - 4 * vstep * Miter, hstep, Prem,
                        @view(μ[hstep * n + 1:hstep * (n + 1)]),
                        @view(σinv[hstep * n + 1:hstep * (n + 1)])
                    )
                end
                if Nrem != 0
                    wait(taskarray[Niter + 1])
                    taskarray[Niter + 1] = @_spawn _ftn!(
                        @view(C[hstep * Niter + 1:end, pstep * Piter + 1:end]),
                        @view(A[vstep * Miter + 1:end, hstep * Niter + 1:end]),
                        @view(B[4 * vstep * Miter + 1:end, pstep * Piter + 1:end]),
                        size(B, 1) - 4 * vstep * Miter, Nrem, Prem,
                        @view(μ[hstep * Niter + 1:end]),
                        @view(σinv[hstep * Niter + 1:end])
                    )
                end
            end
        end
    end
end

for (_ftn!, _ftn_rem!, expr) in [
        (:_snparray_ax_additive!, :_snparray_ax_additive_rem!, 
            :(((Aij >= 2) + (Aij == 3)) * v[j])),
        (:_snparray_ax_dominant!, :_snparray_ax_dominant_rem!, 
            :((Aij >= 2)  * v[j])),
        (:_snparray_ax_recessive!, :_snparray_ax_recessive_rem!, 
            :((Aij == 3) * v[j])),
        (:_snparray_ax_additive_meanimpute!, :_snparray_ax_additive_meanimpute_rem!, 
            :(((Aij >= 2) * 1.0 + (Aij == 3) * 1.0 + (Aij == 1) * μ[j]) * v[j])),
        (:_snparray_ax_dominant_meanimpute!, :_snparray_ax_dominant_meanimpute_rem!, 
            :((Aij >= 2) * v[j] + (Aij == 1) * μ[j] * v[j])),
        (:_snparray_ax_recessive_meanimpute!, :_snparray_ax_recessive_meanimpute_rem!, 
            :((Aij == 3) * v[j] + (Aij == 1) * μ[j] * v[j]))
    ]
    @eval begin
        function ($_ftn_rem!)(out, s, v, μ)
            maxp = length(out)
            @avx for j in eachindex(v)
                block = s[1, j]
                for p in 1:maxp
                    Aij = (block >> (2 * (p - 1))) & 3
                    out[p] += $expr
                end
            end
        end

        function ($_ftn!)(out, s, v, rows, cols, μ)
            k = rows >> 2
            rem = rows & 3

            if k ≥ 1 # avoid avxing over empty collection
                @avx for j ∈ 1:cols
                    for l in 1:k
                        block = s[l, j]

                        for p in 1:4
                            Aij = (block >> (2 * (p - 1))) & 3
                            out[4 * (l - 1) + p] += $expr
                        end

                    end
                end
            end
            if rem != 0
                ($_ftn_rem!)(@view(out[4k + 1:end]), @view(s[k + 1:k + 1, :]), v, μ)
            end
            nothing
        end
    end
end

for (_ftn!, _ftn_rem!, expr) in [
        (:_snparray_atx_additive!, :_snparray_atx_additive_rem!, 
            :(((Aij >= 2) + (Aij == 3)) * v[4 * (l - 1) + p])),
        (:_snparray_atx_dominant!, :_snparray_atx_dominant_rem!, 
            :((Aij >= 2)  * v[4 * (l - 1) + p])),
        (:_snparray_atx_recessive!, :_snparray_atx_recessive_rem!, 
            :((Aij == 3) * v[4 * (l - 1) + p])),
        (:_snparray_atx_additive_meanimpute!, :_snparray_atx_additive_meanimpute_rem!, 
            :(((Aij >= 2) * 1.0 + (Aij == 3) * 1.0 + (Aij == 1) * μ[i]) *  v[4 * (l - 1) + p])),
        (:_snparray_atx_dominant_meanimpute!, :_snparray_atx_dominant_meanimpute_rem!, 
            :((Aij >= 2) * v[4 * (l - 1) + p] + (Aij == 1) * μ[i] * v[4 * (l - 1) + p])),
        (:_snparray_atx_recessive_meanimpute!, :_snparray_atx_recessive_meanimpute_rem!, 
            :((Aij == 3) * v[4 * (l - 1) + p] + (Aij == 1) * μ[i] * v[4 * (l - 1) + p]))
    ]
    @eval begin
        function $(_ftn_rem!)(out, s, v, μ)
            maxp = length(v)
            l = 1
            @avx for i in eachindex(out)
                block = s[1, i]
                for p in 1:maxp
                    Aij = (block >> (2 * (p - 1))) & 3
                    out[i] += $expr
                end
            end
        end

        function $(_ftn!)(out, s, v, rows, cols, μ)
            k = rows >> 2
            rem = rows & 3

            if k ≥ 1 # avoid avxing over empty collection
                @avx for i ∈ 1:cols
                    for l in 1:k
                        block = s[l, i]
                        
                        for p in 1:4
                            Aij = (block >> (2 * (p - 1))) & 3
                            out[i] += $expr
                        end
                    end
                end
            end
            if rem != 0
                $(_ftn_rem!)(out, @view(s[k + 1:k + 1, :]), @view(v[4k + 1:end]), μ)
            end
            nothing
        end
    end
end

for (_ftn!, _ftn_rem!, expr) in [
        (:_snparray_AX_additive!, :_snparray_AX_additive_rem!, 
            :(((Aij >= 2) + (Aij == 3) - μ[j]) * σinv[j] * V[j, c])),
        (:_snparray_AX_dominant!, :_snparray_AX_dominant_rem!, 
            :(((Aij >= 2) - μ[j]) * σinv[j] * V[j, c])),
        (:_snparray_AX_recessive!, :_snparray_AX_recessive_rem!, 
            :(((Aij == 3) - μ[j])* σinv[j] * V[j, c])),
        (:_snparray_AX_additive_meanimpute!, :_snparray_AX_additive_meanimpute_rem!, 
            :(((((Aij >= 2) + (Aij == 3) - μ[j]) + (Aij == 1) * μ[j]) * σinv[j] * V[j, c]))),
        (:_snparray_AX_dominant_meanimpute!, :_snparray_AX_dominant_meanimpute_rem!, 
            :(((Aij >= 2) - μ[j]) * σinv[j] * V[j, c] + (Aij == 1) * μ[j] * σinv[j] * V[j, c])),
            (:_snparray_AX_recessive_meanimpute!, :_snparray_AX_recessive_meanimpute_rem!, 
            :(((Aij == 3) - μ[j]) * V[j, c] * σinv[j] + (Aij == 1) * μ[j] * σinv[j] * V[j, c]))
    ]
    @eval begin
        function ($_ftn_rem!)(out, s, V, μ, σinv, c)
            maxp = size(out, 1)
            @avx for j in 1:size(V, 1)
                block = s[1, j]
                for p in 1:maxp
                    Aij = (block >> (2 * (p - 1))) & 3
                    out[p, c] += $expr
                end
            end
        end

        function ($_ftn!)(out, s, V, srows, scols, Vcols, μ, σinv)
            k = srows >> 2 # fast div(srows, 4)
            rem = srows & 3 # fast rem(srows, 4)

            # compute out[i, c] = s[i, j] * V[j, c] for j in 1:n
            if k ≥ 1 # avoid avxing over empty collection
                @avx for c in 1:Vcols
                    for j in 1:scols
                        for l in 1:k
                            block = s[l, j]
                            for p in 1:4
                                Aij = (block >> (2 * (p - 1))) & 3
                                out[4*(l - 1) + p, c] += $expr
                            end
                        end
                    end
                end
            end
            if rem != 0
                for c in 1:Vcols
                    ($_ftn_rem!)(@view(out[4k + 1:end, :]), @view(s[k + 1:k + 1, :]),
                        V, @view(μ[4k+1:end]), @view(σinv[4k+1:end]), c)
                end
            end
            nothing
        end
    end
end

for (_ftn!, _ftn_rem!, expr) in [
    (:_snparray_AtX_additive!, :_snparray_AtX_additive_rem!, 
        :(((Aij >= 2) + (Aij == 3) - μ[i]) * σinv[i] * V[4 * (l - 1) + p, c])),
    (:_snparray_AtX_dominant!, :_snparray_AtX_dominant_rem!, 
        :(((Aij >= 2) - μ[i]) * σinv[i] * V[4 * (l - 1) + p, c])),
    (:_snparray_AtX_recessive!, :_snparray_AtX_recessive_rem!, 
        :(((Aij == 3) - μ[i]) * σinv[i] * V[4 * (l - 1) + p, c])),
    (:_snparray_AtX_additive_meanimpute!, :_snparray_AtX_additive_meanimpute_rem!, 
        :(((((Aij >= 2) + (Aij == 3) - μ[i]) + (Aij == 1) * μ[i]) * σinv[i] * V[4 * (l - 1) + p, c]))),
    (:_snparray_AtX_dominant_meanimpute!, :_snparray_AtX_dominant_meanimpute_rem!, 
        :(((Aij >= 2) - μ[i]) * σinv[i] * V[4 * (l - 1) + p, c] + (Aij == 1) * μ[i] * σinv[i] * V[4 * (l - 1) + p, c])),
    (:_snparray_AtX_recessive_meanimpute!, :_snparray_AtX_recessive_meanimpute_rem!, 
        :(((Aij == 3) - μ[i]) * σinv[i] * V[4 * (l - 1) + p, c] + (Aij == 1) * μ[i] * σinv[i] * V[4 * (l - 1) + p, c]))
]
    @eval begin
        function ($_ftn_rem!)(out, s, V, μ, σinv, c)
            maxp = size(V, 1)
            l = 1
            @avx for i in 1:size(out, 1)
                block = s[1, i]
                for p in 1:maxp
                    Aij = (block >> (2 * (p - 1))) & 3
                    out[i, c] += $expr
                end
            end
        end

        function ($_ftn!)(out, s, V, srows, scols, Vcols, μ, σinv)
            k = srows >> 2 # fast div(srows, 4)
            rem = srows & 3 # fast rem(srows, 4)

            if k ≥ 1 # avoid avxing over empty collection
                @avx for c in 1:Vcols
                    for i in 1:scols
                        for l in 1:k
                            block = s[l, i]
                            for p in 1:4
                                Aij = (block >> (2 * (p - 1))) & 3
                                out[i, c] += $expr
                            end
                        end
                    end
                end
            end
            if rem != 0
                for c in 1:Vcols
                    ($_ftn_rem!)(out, @view(s[k + 1:k + 1, :]),
                        @view(V[4k + 1:end, :]), μ, σinv, c)
                end
            end
            nothing
        end
    end
end

"""
    Base.copyto!(v, s)

Copy SnpLinAlg `s` to numeric vector or matrix `v`. If `s` is centered/scaled,
`v` will be centered/scaled using precomputed column mean `s.μ` and inverse std 
`s.σinv`.
"""
function Base.copyto!(
    v::AbstractVecOrMat{T}, 
    s::AbstractSnpLinAlg
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

Convert a AbstractSnpLinAlg `s` to a numeric vector or matrix of same shape as `s`.
If `s` is centered/scaled, `v` will be centered/scaled using precomputed column
mean `s.μ` and inverse std `s.σinv`.

# Arguments
- `t::Type{AbstractVecOrMat{T}}`: Vector or matrix type.
"""
Base.convert(::Type{T}, s::AbstractSnpLinAlg) where T <: Array = T(s)
Array{T,N}(s::AbstractSnpLinAlg) where {T,N} = 
    copyto!(Array{T,N}(undef, size(s)), s)
