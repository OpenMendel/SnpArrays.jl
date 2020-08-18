# Hardy-Weinberg: A1A1: (1-p2)^2, A1A2: 2p2(1-p2), A2A2: p2^2
# ADDITIVE_MODEL=1 (A1A1: 0, A1A2: 1, A2A2: 2)
@inline function convert(::Type{T}, a2freq::AbstractFloat, ::Val{1}) where T <: AbstractFloat
    T(2a2freq)
end
# DOMINANT_MODEL=2 (A1A1: 0, A1A2: 1, A2A2: 1)
@inline function convert(::Type{T}, a2freq::AbstractFloat, ::Val{2}) where T <: AbstractFloat
    T(a2freq * (2 - a2freq))
end
# RECESSIVE_MODEL=3 (A1A1: 0, A1A2: 0, A2A2: 1)
@inline function convert(::Type{T}, a2freq::AbstractFloat, ::Val{3}) where T <: AbstractFloat
    T(abs2(a2freq))
end

"""
    Base.copyto!(v, s, P, Q; model=ADDITIVE_MODEL, center=false, scale=false)

Copy SnpArray `s` to a numeric vector or matrix `v`. Impute missing genotypes 
according to ADMIXTURE estimates `P` (population A2 allele frequencies of SNPs) 
and `Q` (population fractions of individuals). Columns of `P` should match columns of 
`s`. Columns of `Q` should match rows of `s`. Note if the inferred minor allele 
frequency from `P` and `Q` is ≤0.01, that genotype will be imputed according to 
the inferred minor allele frequency.

# Positional arguments
- `v`: output vector or array. 
- `s`: SnpArray. 
- `P`: `K x S` array of A2 popluation allele frequencies.
- `Q`: `K x N` array of individual population fractions.

# Keyword arguments
- `model::Union{Val{1}, Val{2}, Val{3}}=ADDITIVE_MODEL`: `ADDITIVE_MODEL` (default), `DOMINANT_MODEL`, or `RECESSIVE_MODEL`.  
- `center::Bool=false`: center each genotype dosage by its mean.
- `scale::Bool=false`: scale each genotype dosage by its standard deviation.
"""
function Base.copyto!(
    v      :: AbstractVecOrMat{T1}, 
    s      :: AbstractSnpArray,
    P      :: AbstractMatrix{T2},
    Q      :: AbstractMatrix{T2};
    model  :: Union{Val{1}, Val{2}, Val{3}} = ADDITIVE_MODEL,
    center :: Bool = false,
    scale  :: Bool = false
    ) where {T1 <: AbstractFloat, T2 <: AbstractFloat}
    m, n, K = size(s, 1), size(s, 2), size(P, 1)
    # center, scale, impute
    @inbounds for j in 1:n, i in 1:m
        vij = SnpArrays.convert(T1, s[i, j], model)
        a2freq = T2(0)
        for k in 1:K
            a2freq += Q[k, i] * P[k, j]
        end
        μij = SnpArrays.convert(T1, a2freq, model)
        σij = model == ADDITIVE_MODEL ? sqrt(μij * (1 - μij / 2)) : sqrt(μij * (1 - μij))
        v[i, j] = (isnan(vij) || a2freq ≤ 0.01 || a2freq ≥ 0.99) ? μij : vij
        center && (v[i, j] -= μij)
        scale  && (v[i, j] /= σij)
    end
    v
end

function Base.convert(
    ::Type{T},
    s::AbstractSnpArray,
    P::AbstractMatrix,
    Q::AbstractMatrix;
    kwargs...) where T <: Array
    T(s, P, Q; kwargs...)
end
Array{T,N}(s::AbstractSnpArray, P::AbstractMatrix, Q::AbstractMatrix; kwargs...) where {T,N} = 
copyto!(Array{T,N}(undef, size(s)), s, P, Q; kwargs...)

function grm_admixture(
    s::AbstractSnpArray,
    P::AbstractMatrix,
    Q::AbstractMatrix,
    ::Type{T} = Float64
    ) where T <: AbstractFloat
    m, n = size(s)
    # convert genotype
    tic = time()
    G = Mmap.mmap(Matrix{T}, m, n) # slightly faster than G = Matrix{T}(undef, m, n)
    Base.copyto!(G, s, P, Q, model=ADDITIVE_MODEL, center=true, scale=true)
    @printf("convert genotype: %.2f seconds\n", time() - tic)
    # GG'
    tic = time()
    Φ = G * transpose(G) # m x m
    @printf("Φ = GG': %.2f seconds\n", time() - tic)
    # convert G to {0,1} matrix
    tic = time()
    @inbounds for (idx, g) in enumerate(G)
        G[idx] = ifelse(iszero(g), T(0), T(1))
        # iszero(g) || (G[idx] = T(1))
    end
    @printf("convert G to {0,1} matrix: %.2f seconds\n", time() - tic)
    # Sij = # SNPs observed in both individuals i and j
    tic = time()
    S = G * transpose(G)
    @printf("S = GG': %.2f seconds\n", time() - tic)
    # Φ /= 2S
    @inbounds for j in 1:m, i in 1:j
        Φ[i, j] /= 2S[i, j]
    end
    copytri!(Φ, 'U')
end # function grm_admixture
