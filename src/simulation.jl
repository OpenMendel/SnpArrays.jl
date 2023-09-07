# Vector of all possible genotypes across one row of compressed SNP data
const genotypes = [0b00_00_00_00;
                   0b00_00_00_10;
                   0b00_00_00_11;
                   0b00_00_10_00;
                   0b00_00_10_10;
                   0b00_00_10_11;
                   0b00_00_11_00;
                   0b00_00_11_10;
                   0b00_00_11_11;
                   0b00_10_00_00;
                   0b00_10_00_10;
                   0b00_10_00_11;
                   0b00_10_10_00;
                   0b00_10_10_10;
                   0b00_10_10_11;
                   0b00_10_11_00;
                   0b00_10_11_10;
                   0b00_10_11_11;
                   0b00_11_00_00;
                   0b00_11_00_10;
                   0b00_11_00_11;
                   0b00_11_10_00;
                   0b00_11_10_10;
                   0b00_11_10_11;
                   0b00_11_11_00;
                   0b00_11_11_10;
                   0b00_11_11_11;
                   0b10_00_00_00;
                   0b10_00_00_10;
                   0b10_00_00_11;
                   0b10_00_10_00;
                   0b10_00_10_10;
                   0b10_00_10_11;
                   0b10_00_11_00;
                   0b10_00_11_10;
                   0b10_00_11_11;
                   0b10_10_00_00;
                   0b10_10_00_10;
                   0b10_10_00_11;
                   0b10_10_10_00;
                   0b10_10_10_10;
                   0b10_10_10_11;
                   0b10_10_11_00;
                   0b10_10_11_10;
                   0b10_10_11_11;
                   0b10_11_00_00;
                   0b10_11_00_10;
                   0b10_11_00_11;
                   0b10_11_10_00;
                   0b10_11_10_10;
                   0b10_11_10_11;
                   0b10_11_11_00;
                   0b10_11_11_10;
                   0b10_11_11_11;
                   0b11_00_00_00;
                   0b11_00_00_10;
                   0b11_00_00_11;
                   0b11_00_10_00;
                   0b11_00_10_10;
                   0b11_00_10_11;
                   0b11_00_11_00;
                   0b11_00_11_10;
                   0b11_00_11_11;
                   0b11_10_00_00;
                   0b11_10_00_10;
                   0b11_10_00_11;
                   0b11_10_10_00;
                   0b11_10_10_10;
                   0b11_10_10_11;
                   0b11_10_11_00;
                   0b11_10_11_10;
                   0b11_10_11_11;
                   0b11_11_00_00;
                   0b11_11_00_10;
                   0b11_11_00_11;
                   0b11_11_10_00;
                   0b11_11_10_10;
                   0b11_11_10_11;
                   0b11_11_11_00;
                   0b11_11_11_10;
                   0b11_11_11_11
                   ];

missing_mask    = [0b00_00_00_00;
                   0b00_00_00_11;
                   0b00_00_11_00;
                   0b00_00_11_11;
                   0b00_11_00_00;
                   0b00_11_00_11;
                   0b00_11_11_00;
                   0b00_11_11_11;
                   0b11_00_00_00;
                   0b11_00_00_11;
                   0b11_00_11_00;
                   0b11_00_11_11;
                   0b11_11_00_00;
                   0b11_11_00_11;
                   0b11_11_11_00;
                   0b11_11_11_11;
                   ];

"""
    compressed_snp_prob!(compressed_prob, prob)

Populate `compressed_prob` with probabilities of each possible observed value in 
a compressed SNP format, assuming non-missing genotype data
"""
@inline function compressed_snp_prob!(
    compressed_prob::AbstractVector{T}, 
    prob::T
    ) where {T<:AbstractFloat}
    # Tuple is non-allocating since its immutable
    probs = (abs2(1 - prob), 2 * prob * (1 - prob), abs2(prob))
    @inbounds for l in 1:3
        for k in 1:3
            for j in 1:3
                for i in 1:3
                    idx = (l - 1) * 27 + (k - 1) * 9 + (j - 1) * 3 + i
                    compressed_prob[idx] = probs[l] * probs[k] * probs[j] * probs[i] 
                end
            end
        end
    end
    return compressed_prob
end

@inline function compressed_missing_prob!(
    compressed_prob::AbstractVector{T}, 
    prob::T
    ) where {T<:AbstractFloat}
    probs = (1 - prob, prob)
    @inbounds for l in 1:2
        for k in 1:2
            for j in 1:2
                for i in 1:2
                    idx = (l - 1) * 8 + (k - 1) * 4 + (j - 1) * 2 + i
                    compressed_prob[idx] = probs[l] * probs[k] * probs[j] * probs[i] 
                end
            end
        end
    end
    return compressed_prob
end

struct compressedSNPWeights{S<:AbstractFloat, T<:AbstractFloat, V<:AbstractVector{T}} <: StatsBase.AbstractWeights{S, T, V}
    values::V
    sum::S
    function compressedSNPWeights{S, T, V}(values, sum) where {S<:AbstractFloat, T<:AbstractFloat, V<:AbstractVector{T}}
        isfinite(sum) || throw(ArgumentError("weights cannot contain Inf or NaN values"))
        return new{S, T, V}(values, sum)
    end
    compressedSNPWeights(values::AbstractVector{T}, sum::S) where {S<:AbstractFloat, T<:AbstractFloat} = compressedSNPWeights{S, T, typeof(values)}(values, sum)
end

function compressedSNPWeights(values::AbstractVector{T}) where {T<:AbstractFloat}
    return compressedSNPWeights(values, one(T))
end

struct compressedMissingWeights{S<:AbstractFloat, T<:AbstractFloat, V<:AbstractVector{T}} <: StatsBase.AbstractWeights{S, T, V}
    values::V
    sum::S
    function compressedMissingWeights{S, T, V}(values, sum) where {S<:AbstractFloat, T<:AbstractFloat, V<:AbstractVector{T}}
        isfinite(sum) || throw(ArgumentError("weights cannot contain Inf or NaN values"))
        return new{S, T, V}(values, sum)
    end
    compressedMissingWeights(values::AbstractVector{T}, sum::S) where {S<:AbstractFloat, T<:AbstractFloat} = compressedMissingWeights{S, T, typeof(values)}(values, sum)
end

function compressedMissingWeights(values::AbstractVector{T}) where {T<:AbstractFloat}
    return compressedMissingWeights(values, one(T))
end

@inline function _to_missing!(g::AbstractVector{UInt8}, mask::AbstractVector{UInt8})
    @inbounds @simd for i in eachindex(g, mask)
        g[i] = (g[i] & ~mask[i]) | (mask[i] & 0x55)
    end
    return g
end

function _make_alias_table!(
    w      :: AbstractVector{T}, 
    wsum   :: T,
    a      :: AbstractVector{Float64},
    alias  :: AbstractVector{Int};
    larges :: Vector{Int} = Vector{Int}(undef, length(w)),
    smalls :: Vector{Int} = Vector{Int}(undef, length(w))
    ) where {T<:AbstractFloat}
    # Arguments:
    #
    #   w [in]:         input weights
    #   wsum [in]:      pre-computed sum(w)
    #
    #   a [out]:        acceptance probabilities
    #   alias [out]:    alias table
    #
    # Note: a and w can be the same array, then that array will be
    #       overwritten inplace by acceptance probabilities
    #
    # Returns nothing

    n = length(w)
    length(a) == length(alias) == n ||
    throw(DimensionMismatch("Inconsistent array lengths."))

    # For our use-case, wsum is always equal to 1.0
    ac = n / wsum
    for i in 1:n
        @inbounds a[i] = w[i] * ac
    end

    kl = 0  # actual number of larges
    ks = 0  # actual number of smalls

    for i in 1:n
        @inbounds ai = a[i]
        if ai > 1.0
            larges[kl += 1] = i  # push to larges
        elseif ai < 1.0
            smalls[ks += 1] = i  # push to smalls
        end
    end

    while kl > 0 && ks > 0
        s = smalls[ks]; ks -= 1  # pop from smalls
        l = larges[kl]; kl -= 1  # pop from larges
        @inbounds alias[s] = l
        @inbounds al = a[l] = (a[l] - 1.0) + a[s]
        if al > 1.0
            larges[kl += 1] = l  # push to larges
        else
            smalls[ks += 1] = l  # push to smalls
        end
    end

    # this loop should be redundant, except for rounding
    for i in 1:ks
        @inbounds a[smalls[i]] = 1.0
    end
    return nothing
end

# Assumes alias table has been pre-computed and stored in ap and alias
function _alias_sample!(
    rng   :: Random.AbstractRNG,
    a     :: AbstractArray{T},
    x     :: AbstractArray{T},
    ap    :: AbstractVector{Float64},
    alias :: AbstractVector{Int}
    ) where {T<:Real}

    Base.mightalias(a, x) &&
        throw(ArgumentError("output array x must not share memory with input array a"))
    1 == firstindex(a) == firstindex(x) ||
        throw(ArgumentError("non 1-based arrays are not supported"))

    n = length(a)

    # sampling
    S = Random.Sampler(rng, 1:n)
    for i in eachindex(x)
        j = rand(rng, S)
        x[i] = rand(rng) < ap[j] ? a[j] : a[alias[j]]
    end
    return x
end

function _simulate!(
    rng::Random.AbstractRNG,
    X::AbstractMatrix{UInt8},
    mafs::AbstractVector{T},
    weights::compressedSNPWeights
    ) where {T<:AbstractFloat}
    # pre-allocate vectors for alias table sampling
    larges = Vector{Int}(undef, 81)
    smalls = Vector{Int}(undef, 81)

    ap    = Vector{T}(undef, 81)
    alias = Vector{Int}(undef, 81)

    @inbounds for j in axes(X, 2)
        ρ          = mafs[j]
        # # Store probabilities for 00, 01, and 11
        # storage[1] = abs2(1 - ρ)
        # storage[2] = 2 * ρ * (1 - ρ)
        # storage[3] = abs2(ρ)
        # Map probabilities for one observation to probabilities for 4
        compressed_snp_prob!(weights.values, ρ)
        # Construct alias table takes O(k log(k)) where k = 81
        _make_alias_table!(weights.values, one(T), ap, alias; larges=larges, smalls=smalls)
        # Drawing samples is only O(1) per sample, should be fastest method
        _alias_sample!(rng, genotypes, view(X, :, j), ap, alias)
    end
    return X
end

function _simulate_missing!(
    rng::Random.AbstractRNG,
    X::AbstractMatrix{UInt8},
    missing_rates::AbstractVector{T},
    weights::compressedMissingWeights
    ) where {T<:AbstractFloat}
    larges = Vector{Int}(undef, 16)
    smalls = Vector{Int}(undef, 16)

    ap = Vector{Float64}(undef, 16)
    alias = Vector{Int}(undef, 16)

    missing_idx = Vector{UInt8}(undef, size(X, 1))
    @inbounds for j in axes(X, 2)
        ρ = missing_rates[j]
        # Map missing probability for one observation to probabilities for 4
        compressed_missing_prob!(weights.values, ρ)
        # Construct alias table takes O(k log(k)) where k = 81
        _make_alias_table!(weights.values, one(T), ap, alias; larges=larges, smalls=smalls)
        # Drawing samples is only O(1) per sample, should be fastest method
        # Sample index vector to mask off missing values
        _alias_sample!(rng, missing_mask, missing_idx, ap, alias)
        # Fast operation to assign missing according to indexed values
        _to_missing!(view(X, :, j), missing_idx)
    end
    return X
end

"""
    simulate!([rng], m, n, mafs)

Simulate genotype data to be stored directly in a `SnpArray`, assuming no missingness.

# Arguments
- `m::Integer`: Number of subjects
- `n::Integer`: Number of columns
- `mafs::AbstractVector{T}`: Vector of minor allele frequencies, should be length `n`
"""

function simulate!(
    rng::Random.AbstractRNG,
    m::Integer, 
    n::Integer,
    mafs::AbstractVector{T}
    ) where {T<:AbstractFloat}
    length(mafs) == n || throw(DimensionMismatch("Number of MAFs must equal number of columns"))

    # Get number of expected rows in compressed SNP data
    r = cld(m, 4)

    weights = compressedSNPWeights(Vector{T}(undef, 81))
    # storage = Vector{T}(undef, 3)

    # Allocate storage
    X_compressed = Matrix{UInt8}(undef, r, n)

    # Simulate data
    _simulate!(rng, X_compressed, mafs, weights)

    return SnpArray(X_compressed, zeros(Int, (4, n)), zeros(Int, (4, m)), m)
end

function simulate!(
    m::Integer, 
    n::Integer,
    mafs::AbstractVector{T}
    ) where {T<:AbstractFloat}
    return simulate!(Random.GLOBAL_RNG, m, n, mafs)
end

"""
    simulate!([rng], X, mafs)

Simulate genotype data to be stored directly in a `SnpArray`, assuming no missingness.

# Arguments
- `X::SnpArray`: `m` by `n` SnpArray
- `mafs::AbstractVector{T}`: Vector of minor allele frequencies, should be length `n`
"""

function simulate!(
    rng::Random.AbstractRNG,
    X::SnpArray,
    mafs::AbstractVector{T}
    ) where {T<:AbstractFloat}
    m, n = size(X)
    length(mafs) == n || throw(DimensionMismatch("Number of MAFs must equal number of columns"))

    weights = compressedSNPWeights(Vector{T}(undef, 81))

    # Simulate data
    _simulate!(rng, X.data, mafs, weights)

    return X
end

function simulate!(
    X::SnpArray,
    mafs::AbstractVector{T}
    ) where {T<:AbstractFloat}
    return simulate!(Random.GLOBAL_RNG, X, mafs)
end

"""
    simulate_missing!([rng], X, missing_rates)

In-place modify a `SnpArray` to simulate missing data

# Arguments
- `X::SnpArray`: `m` by `n` SnpArray
- `missing_rates::AbstractVector{T}`: Vector of missingness rates, should be length `n`
"""

function simulate_missing!(
    rng::Random.AbstractRNG,
    X::SnpArray,
    missing_rates::AbstractVector{T}
    ) where {T<:AbstractFloat}
    m, n = size(X)
    length(missing_rates) == n || throw(DimensionMismatch("Number of missing rates must equal number of columns"))

    weights = compressedMissingWeights(Vector{T}(undef, 16))

    # Simulate data
    _simulate_missing!(rng, X.data, missing_rates, weights)
    return X
end

function simulate_missing!(
    X::SnpArray,
    missing_rates::AbstractVector{T}
    ) where {T<:AbstractFloat}
    return simulate_missing!(Random.GLOBAL_RNG, X, missing_rates)
end