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
    if model == ADDITIVE_MODEL
        if center || scale || impute
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
        if center || scale || impute
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
        if center || scale || impute
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
    LinearAlgebra.mul!(out, s::SnpLinAlg, v)

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
    LinearAlgebra.mul!(out, s::SnpLinAlg, v)

In-place matrix-matrix multiplication.
"""
function mul!(
    out::AbstractMatrix{T}, 
    sla::SnpLinAlg{T}, 
    v::AbstractMatrix{T}) where T <: AbstractFloat
    @assert size(out, 1) == size(sla, 1) && size(out, 2) == size(v, 2) && size(sla, 2) == size(v, 1)

    fill!(out, zero(eltype(out)))
    s = sla.s

    _snparray_AX_tile!(out, s.data, w, sla.model, sla.μ, sla.impute)
end

"""
    LinearAlgebra.mul!(out, s::Union{Transpose{T, SnpLinAlg{T}}, Adjoint{T, SnpLinAlg{T}}}, v)

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
    Mrem = M & (vstep - 1) # fast rem(M, 512)
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
                    vstep << 2, hstep, μ
                )
            end
            if Mrem != 0
                wait(taskarray[Miter+1])
                taskarray[Miter+1] = @_spawn _ftn!(
                    @view(c[4 * vstep * Miter + 1:end]), 
                    @view(A[vstep * Miter + 1:end, hstep * n + 1:hstep * (n + 1)]),
                    @view(b[hstep * n + 1:hstep * (n + 1)]),
                    length(c) - 4 * vstep * Miter, hstep, μ
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
                    Nrem, μ
                )
            end
            if Mrem != 0
                wait(taskarray[Miter + 1])
                taskarray[Miter + 1] = @_spawn _ftn!(
                    @view(c[4 * vstep * Miter+1:end]),
                    @view(A[vstep * Miter + 1:end, hstep * Niter + 1:end]),
                    @view(b[hstep * Niter + 1:end]),
                    length(c) - 4 * vstep * Miter,
                    Nrem, μ
                )
            end
        end
    end
end

function _snparray_AX_tile!(C, A, B, model, μ, impute)
    vstep = 1024
    hstep = 1024
    vstep_log2 = 10
    hstep_log2 = 10

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
    
    fill!(C, zero(UInt8))

    M = length(c) >> 2
    N = size(A, 2)
    P = size(C, 2)
    Miter = M >>> vstep_log2 # fast div(M, vstep_log2)
    Mrem = M & (vstep - 1) # fast rem(M, vstep)
    Niter = N >>> hstep_log2
    Nrem = N & (hstep - 1)
    Piter = P >>> hstep_log2
    Prem = P & (hstep - 1)
    taskarray = Array{Any}(undef, Miter + 1)
    fill!(taskarray, nothing)
    @_sync begin
        GC.@preserve C A B for n in 0:Niter - 1
            for m in 0:Miter - 1
                for p in 0:Piter - 1
                    wait(taskarray[m+1])
                    taskarray[m+1] = @_spawn _ftn!(
                        gesp(stridedpointer(C), (4 * vstep * m, pstep * p)),
                        gesp(stridedpointer(A), (vstep * m, hstep * n)),
                        gesp(stridedpointer(B), (hstep * n, pstep * p)),
                        vstep << 2, hstep, pstep, μ
                    )
                end
            end
            if Mrem != 0
                wait(taskarray[Miter+1])
                taskarray[Miter+1] = @_spawn _ftn!(
                    @view(C[4 * vstep * Miter + 1:end]), 
                    @view(A[vstep * Miter + 1:end, hstep * n + 1:hstep * (n + 1)]),
                    @view(B[hstep * n + 1:hstep * (n + 1)]),
                    length(c) - 4 * vstep * Miter, hstep, μ
                )
            end
        end
        if Nrem != 0
            for m in 0:Miter-1
                wait(taskarray[m+1])
                taskarray[m+1] = @_spawn _ftn!(
                    @view(C[4 * vstep * m + 1:4 * vstep * (m + 1)]),
                    @view(A[vstep * m + 1:vstep * (m + 1), hstep * Niter + 1:end]),
                    @view(BLAS[hstep * Niter + 1:end]),
                    vstep << 2, 
                    Nrem, μ
                )
            end
            if Mrem != 0
                wait(taskarray[Miter + 1])
                taskarray[Miter + 1] = @_spawn _ftn!(
                    @view(C[4 * vstep * Miter+1:end]),
                    @view(A[vstep * Miter + 1:end, hstep * Niter + 1:end]),
                    @view(BitSet[hstep * Niter + 1:end]),
                    length(c) - 4 * vstep * Miter,
                    Nrem, μ
                )
            end
        end
    end
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
    Mrem = M & (vstep - 1) # fast rem(M, 512)
    Niter = N >>> hstep_log2
    Nrem = N & (hstep - 1)

    taskarray = Array{Any}(undef, Niter+1)
    fill!(taskarray, nothing)
    @_sync begin
        GC.@preserve c A b for m in 0:Miter - 1
            for n in 0:Niter - 1
                wait(taskarray[n + 1])
                @_spawn _ftn!(
                    gesp(stridedpointer(c), (hstep * n,)),
                    gesp(stridedpointer(A), (vstep * m, hstep * n)),
                    gesp(stridedpointer(b), (4 * vstep * m,)),
                    vstep << 2, hstep, μ
                )
            end
            if Nrem != 0
                wait(taskarray[Niter + 1])
                @_spawn _ftn!(
                    @view(c[hstep * Niter + 1:end]),
                    @view(A[vstep * m + 1:vstep * (m + 1), hstep * Niter + 1:end]),
                    @view(b[4 * vstep * m + 1:4 * vstep * (m + 1)]),
                    vstep << 2, Nrem, μ
                )
            end

        end
        if Mrem != 0
            for n in 0:Niter - 1
                wait(taskarray[n + 1])
                @_spawn _ftn!(
                    @view(c[hstep * n + 1:hstep * (n + 1)]),
                    @view(A[vstep * Miter + 1:end, hstep * n + 1:hstep * (n + 1)]),
                    @view(b[4 * vstep * Miter + 1:end]),
                    length(b) - 4 * vstep * Miter,
                    hstep, μ
                )
            end
            if Nrem != 0
                wait(taskarray[Niter + 1])
                @_spawn _ftn!(
                    @view(c[hstep * Niter + 1:end]),
                    @view(A[vstep * Miter + 1:end, hstep * Niter + 1:end]),
                    @view(b[4 * vstep * Miter + 1:end]),
                    length(b) - 4 * vstep * Miter,
                    Nrem, μ
                )
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
            :(((Aij >= 2) * 1.0 + (Aij == 3) * 1.0 + (Aij == 1) * μ[j]) *  v[j])),
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

            @avx for j ∈ 1:cols
                for l in 1:k
                    block = s[l, j]

                    for p in 1:4
                        Aij = (block >> (2 * (p - 1))) & 3
                        out[4 * (l - 1) + p] += $expr
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

            @avx for i ∈ 1:cols
                for l in 1:k
                    block = s[l, i]
                    
                    for p in 1:4
                        Aij = (block >> (2 * (p - 1))) & 3
                        out[i] += $expr
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
            :(((Aij >= 2) + (Aij == 3)) * V[j, c])),
        (:_snparray_AX_dominant!, :_snparray_AX_dominant_rem!, 
            :((Aij >= 2)  * V[j, c])),
        (:_snparray_AX_recessive!, :_snparray_AX_recessive_rem!, 
            :((Aij == 3) * V[j, c])),
        (:_snparray_AX_additive_meanimpute!, :_snparray_AX_additive_meanimpute_rem!, 
            :(((Aij >= 2) * 1.0 + (Aij == 3) * 1.0 + (Aij == 1) * μ[j]) *  V[j, c])),
        (:_snparray_AX_dominant_meanimpute!, :_snparray_AX_dominant_meanimpute_rem!, 
            :((Aij >= 2) * V[j, c] + (Aij == 1) * μ[j] * V[j, c])),
        (:_snparray_AX_recessive_meanimpute!, :_snparray_AX_recessive_meanimpute_rem!, 
            :((Aij == 3) * V[j, c] + (Aij == 1) * μ[j] * V[j, c]))
    ]
    @eval begin
        function ($_ftn_rem!)(out, s, V, μ)
            maxp = length(out)
            @avx for j in eachindex(V)
                block = s[1, j]
                for p in 1:maxp
                    Aij = (block >> (2 * (p - 1))) & 3
                    out[p] += $expr
                end
            end
        end

        function ($_ftn!)(out, s, V, srows, scols, Vcols, μ)
            k = srows >> 2 # fast div(srows, 4)
            rem = srows & 3 # fast rem(srows, 4)

            # out[i, c] = s[i, j] * V[j, c] for j in 1:scols
            @avx for j in 1:scols
                for l in 1:k
                    block = s[l, j]
                    for c in 1:Vcols
                        for p in 1:4
                            Aij = (block >> (2 * (p - 1))) & 3
                            out[4 * (l - 1) + p, c] += $expr
                        end
                    end

                end
            end
            if rem != 0
                ($_ftn_rem!)(@view(out[4k + 1:end]), @view(s[k + 1:k + 1, :]), V, μ)
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
