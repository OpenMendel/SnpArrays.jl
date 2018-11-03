"""
    grm(s, method=:GRM, minmaf=0.01, colinds=nothing)

Compute empirical kinship matrix from a SnpArray `s`. Missing genotypes are imputed
on the fly by column mean.

# Arguments
- `method::Symbol=:GRM`: `:GRM` (default), `:MoM`, or `:Robust`.
- `minmaf::Real=0.01`: columns with MAF less than `minmaf` are excluded; default 0.01.
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
    # G = Mmap.mmap(Matrix{T}, m, n) # about same speed
    G = Matrix{T}(undef, m, n)
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
