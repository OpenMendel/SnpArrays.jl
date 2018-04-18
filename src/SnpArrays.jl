module SnpArrays

using Compat, LinearMaps

import Base: filter, issymmetric

export estimatesize, filter, grm, _grm, _mom, pca, pca_sp, randgeno,
  SnpArray, SnpData, SnpLike, summarize, writeplink, issym, issymmetric

type SnpArray{N} <: AbstractArray{NTuple{2,Bool}, N}
  A1::BitArray{N}
  A2::BitArray{N}
end

# SnpArray or a view of a SnpArray
@compat SnpLike{N, S<:SnpArray} = Union{SnpArray{N}, SubArray{NTuple{2, Bool}, N, S}}
@compat SnpMatrix = SnpArray{2}
@compat SnpVector = SnpArray{1}

"""
Construct a SnpArray from an array of A1 allele counts {0, 1, 2}.
"""
function SnpArray{T <: Real}(a1count::AbstractArray{T})
    SnpArray(a1count .> one(T), a1count .> zero(T))
end

"""
Create a SnpArray with all A1 alleles.
"""
function SnpArray(dims...)
  SnpArray(falses(dims), falses(dims))
end #

"""
    SnpArray(plinkFile, [people], [snps])

Construct a SnpArray from Plink binary file `plinkFile.bed`.
"""
function SnpArray(
  plinkFile::AbstractString;
  people::Int = countlines(plinkFile * ".fam"),
  snps::Int = countlines(plinkFile * ".bim")
  )

  # dimensions
  (people > 0 && snps > 0) || throw(ArgumentError("people and snps have to be positive integer"))
  # read binary genotype data from bed file
  plinkBedfile = contains(plinkFile, ".bed")? plinkFile : string(plinkFile, ".bed")
  fid = open(plinkBedfile, "r")
  # decide BED file version and snp/individual-major
  bedheader = read(fid, UInt8, 3)
  # PLINK coding (genotype->bits): A1/A1->00, A1/A2->01, A2/A2->11, missing->10
  # http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml
  # v1.0 BED file: first two magic byte, 3rd byte indicate SNP- or individual-major
  if bits(bedheader[1]) == "01101100" && bits(bedheader[2]) == "00011011" &&
    bits(bedheader[3]) == "00000001"
    info("v1.0 BED file detected")
    snp_major = true
    seek(fid, 3)
  elseif bits(bedheader[1]) == "01101100" && bits(bedheader[2]) == "00011011" &&
    bits(bedheader[3]) == "00000000"
    info("v1.0 BED file detected")
    snp_major = false
    seek(fid, 3)
  # v0.99 BED file: no 2-byte magic number, 1st byte indicate SNP- or individual-major
  elseif bits(bedheader[1]) == "00000001"
    info("v0.99 BED file detected")
    snp_major = true
    seek(fid, 1)
  elseif bits(bedheader[1]) == "00000000"
    snp_major = false
    seek(fid, 1)
  # Prior to v0.99: no 2-byte magic number, no SNP-major/individual-major
  # identifier, always individual-major
  else
    info("prior to v0.99 BED file detected")
    snp_major = false
    seek(fid, 0)
  end
  if snp_major
    plinkbits = Mmap.mmap(fid, BitArray, (2, 4ceil(Int, 0.25people), snps))
    A1 = plinkbits[1, 1:people, :]
    A2 = plinkbits[2, 1:people, :]
  else # individual-major
    plinkbits = Mmap.mmap(fid, BitArray, (2, 4ceil(Int, 0.25snps), people))
    A1 = plinkbits[1, 1:snps, :]'
    A2 = plinkbits[2, 1:snps, :]'
  end
  close(fid)
  SnpArray(A1, A2)
end # function SnpArray

#---------------------------------------------------------------------------# methods
# Julia docs on methods required for AbstractArray:
# http://docs.julialang.org/en/release-0.4/manual/interfaces/#man-interfaces-abstractarray

Base.size(A::SnpArray)                 = size(A.A1)
Base.size(A::SnpArray, d::Int)         = size(A.A1, d)
Base.ndims(A::SnpArray)                = ndims(A.A1)
Base.endof(A::SnpArray)                = length(A)
Base.eltype(A::SnpArray)               = NTuple{2, Bool}
@compat Base.IndexStyle(::Type{SnpArray})      = Base.IndexLinear()

@inline function Base.getindex(A::SnpArray, i::Int)
  (getindex(A.A1, i), getindex(A.A2, i))
end

@inline function Base.getindex(A::SnpArray, i::Int, j::Int)
  (getindex(A.A1, i, j), getindex(A.A2, i, j))
end

@inline function Base.setindex!(A::SnpArray, v::NTuple{2, Bool}, i::Int)
  setindex!(A.A1, v[1], i), setindex!(A.A2, v[2], i)
end

@inline function Base.setindex!(A::SnpArray, v::NTuple{2, Bool}, i::Int, j::Int)
  setindex!(A.A1, v[1], i, j), setindex!(A.A2, v[2], i, j)
end

@inline function Base.setindex!(A::SnpArray, v::Real, i::Int)
  # real number v is interpreted as A1 allele count
  if isnan(v)
    setindex!(A, (true, false), i)
  else
    setindex!(A, (v > one(eltype(v)), v > zero(eltype(v))), i)
  end
end

@inline function Base.setindex!(A::SnpArray, v::Real, i::Int, j::Int)
  if isnan(v)
    setindex!(A, (true, false), i, j)
  else
    setindex!(A, (v > one(eltype(v)), v > zero(eltype(v))), i, j)
  end
end

function Base.similar(A::SnpArray, dims::Dims)
  SnpArray(BitArray(dims), BitArray(dims))
end # function Base.similar
#
function Base.similar(A::SnpArray, ::Type{NTuple{2, Bool}}, dims::Dims)
  similar(A, dims)
end # function Base.similar

#---------------------------------------------------------------------------#
# Code for missing genotypes
#---------------------------------------------------------------------------#

# missing code is 10 = (true, false)
@inline Base.isnan(a1::Bool, a2::Bool) = a1 & !a2
@inline Base.isnan(a::Tuple{Bool, Bool}) = Base.isnan(a[1], a[2])
function Base.isnan{N}(A::SnpLike{N})
  b = BitArray(size(A))
  @inbounds @simd for i in eachindex(A)
    b[i] = Base.isnan(A[i])
  end
  return b
end # function Base.isnan

#---------------------------------------------------------------------------#
# Convert and copy
#---------------------------------------------------------------------------#

"""
    estimatesize(n, p, t[, maf])

Estimate the memory usage (in bytes) for storing SNP data with `n` individuals,
`p` SNPs, and average minor allele frequency `maf`. `maf` is used only when
`t <: AbstractSparseMatrix`.
"""
function estimatesize{T <: AbstractFloat}(n::Int, p::Int,
  t::Type, maf::T = 0.25)
  if t <: DenseMatrix
    storage = convert(Float64, sizeof(eltype(t)) * n * p)
  elseif t <: AbstractSparseMatrix
    storage = convert(Float64,
      (sizeof(t.parameters[1]) + sizeof(t.parameters[2])) *
      n * p * maf * (2.0 - maf) + sizeof(t.parameters[2]) * (p + 1))
  end
  storage
end # function estimatesize


"""
Create a SnpArray with all A1 alleles.
"""
function Base.convert(t::Type{SnpArray}, dims...)
  SnpArray(falses(dims), falses(dims))
end # function Base.convert

"""
Convert a two-bit genotype to a real number (minor allele count) according to
specified SNP model. Missing genotype is converted to `NaN`. `minor_allele=true`
indicates `A1` is the minor allele; `minor_allele=false` indicates `A2` is the
minor allele.
"""
function Base.convert{T <: Real}(t::Type{T}, a::NTuple{2, Bool},
  minor_allele::Bool, model::Symbol = :additive)
  if isnan(a)
    b = convert(T, NaN)
  else
    if minor_allele # A1 is the minor allele
      if model == :additive
        b = convert(T, !a[1] + !a[2])
      elseif model == :dominant
        b = convert(T, !a[1] | !a[2])
      elseif model == :recessive
        b = convert(T, !a[1] & !a[2])
      end
    else # A2 is the minor allele
      if model == :additive
        b = convert(T, a[1] + a[2])
      elseif model == :dominant
        b = convert(T, a[1] | a[2])
      elseif model == :recessive
        b = convert(T, a[1] & a[2])
      end
    end
  end
  return b
end # function Base.convert

"""
Convert a SNP matrix to a numeric matrix of minor allele counts according to
specified SNP model. If `impute == true`, missing entries are imputed according
to (column) minor allele frequencies.
"""
function Base.convert{T <: Real, N}(t::Type{Array{T, N}}, A::SnpLike{N};
  model::Symbol = :additive, impute::Bool = false, center::Bool = false,
  scale::Bool = false)
  B = similar(A, T)
  copy!(B, A; model = model, impute = impute, center = center, scale = scale)
end # function Base.convert

function Base.copy!{T <: Real, N}(B::AbstractArray{T, N}, A::SnpLike{N};
  model::Symbol = :additive, impute::Bool = false, center::Bool = false,
  scale::Bool = false)
  size(B) == size(A) || throw(ArgumentError("Dimensions do not match"))
  if ndims(A) == 1
    m, n = length(A), 1
  elseif ndims(A) == 2
    m, n = size(A)
  end
  # convert column by column
  @inbounds for j in 1:n
    # first pass: find minor allele and its frequency
    if ndims(A) == 1
      maf, minor_allele, = summarize(A)
    else ndims(A) == 2
      maf, minor_allele, = summarize(view(A, :, j))
    end
    # second pass: impute, convert, center, scale
    ct = convert(T, 2.0maf)
    wt = convert(T, maf == 0.0 ? 1.0 : 1.0 / √(2.0maf * (1.0 - maf)))
    @simd for i in 1:m
      (a1, a2) = A[i, j]
      # impute if asked
      if isnan(a1, a2) && impute
        a1, a2 = randgeno(maf, minor_allele)
      end
      B[i, j] = convert(T, (a1, a2), minor_allele, model)
      if center; B[i, j] -= ct; end
      if scale; B[i, j] *= wt; end
    end
  end
  return B
end # function Base.copy!

"""
Convert a SNP matrix to a numeric matrix of minor allele counts according to
specified SNP model. If `impute == true`, missing entries are imputed according
to (column) minor allele frequencies.
"""
function Base.convert{T <: Real, TI <: Integer}(t::Type{SparseMatrixCSC{T, TI}},
  A::SnpLike{2}; model::Symbol = :additive, impute::Bool = false)
  m, n = size(A)
  # prepare sparese matrix data structure
  rowval = TI[]
  nzval = T[]
  colptr = zeros(TI, n + 1)
  # convert column by column
  zeroT = convert(T, 0.0)
  @inbounds for j in 1:n
    colptr[j] = convert(TI, length(nzval) + 1)
    # first pass: find minor allele and its frequency
    maf, minor_allele, = summarize(view(A, :, j))
    # second pass: impute, convert
    for i in 1:m
      (a1, a2) = A[i, j]
      # impute if asked
      if isnan(a1, a2)
        if impute
          a1, a2 = randgeno(maf, minor_allele)
        else
          continue # nan and no impute -> done
        end
      end
      v = convert(T, (a1, a2), minor_allele, model)
      if v ≠ zeroT
        push!(rowval, convert(TI, i)), push!(nzval, v)
      end
    end # i
  end # j
  colptr[n + 1] = length(nzval) + 1
  return SparseMatrixCSC(m, n, colptr, rowval, nzval)
end # function Base.convert

#---------------------------------------------------------------------------#
# Random genotype generation
#---------------------------------------------------------------------------#

"""
Generate a genotype according to a1 allele frequency.
"""
function randgeno{T <: AbstractFloat}(a1freq::T)
  b1 = rand() > a1freq
  b2 = rand() > a1freq
  # make sure not conflict with missing code 10
  if isnan(b1, b2)
    (b1, b2) = (false, true)
  end
  return b1, b2
end # function readgeno

"""
Generate a genotype according to minor allele frequency. `minor_allele` indicates
the minor allele is A1 (`true`) or A2 (`false`).
"""
@inline function randgeno{T <: AbstractFloat}(maf::T, minor_allele::Bool)
  minor_allele ? randgeno(maf) : randgeno(1.0 - maf)
end # function readgeno

"""
Generate a SnpVector according to minor allele frequency. `minor_allele`
indicates the minor allele is A1 (`true`) or A2 (`false`).
"""
function randgeno{T <: AbstractFloat}(n::Int, maf::T, minor_allele::Bool)
  s = SnpArray(n)
  @inbounds @simd for i in 1:n
    s[i] = randgeno(maf, minor_allele)
  end
  return s
end # function readgeno

"""
Generate a SnpMatrix according to minor allele frequencies `maf`. `minor_allele`
vector indicates the minor alleles are A1 (`true`) or A2 (`false`).
"""
function randgeno{T <: AbstractFloat}(m::Int, n::Int, maf::Vector{T},
  minor_allele::BitVector)
  length(maf) == n && throw(ArgumentError("length of maf should be n"))
  length(minor_allele) == n && throw(ArgumentError("length of minor_allele should be n"))
  s = SnpArray(m, n)
  @inbounds @simd for j in 1:n
    for i in 1:m
      s[i, j] = randgeno(maf[j], minor_allele[j])
    end
  end
  return s
end # function readgeno

function randgeno{T <: AbstractFloat}(n::Tuple{Int,Int}, maf::Vector{T},
  minor_allele::BitArray{1})
  randgeno(n..., maf, minor_allele)
end # function readgeno

#---------------------------------------------------------------------------#
# Summary statistics
#---------------------------------------------------------------------------#

"""
Compute summary statistics of a SnpMatrix.

# Output:
* `maf` - minor allele frequency of each SNP
* `minor_allele` - indicate the minor allele is A1 (`true`) or A2 (`false`)
* `missings_by_person` - number of missing genotypes for each person
* `missings_by_snp` - number of missing genotypes for each SNP
"""
function summarize(A::SnpLike{2})
  m, n = size(A)
  maf = zeros(Float64, n)             # minor allele frequencies for each column
  minor_allele = trues(n)             # true->A1 is the minor allele
  missings_by_snp = zeros(Int, n)     # no. missing genotypes for each row
  missings_by_person = zeros(Int, m)  # no. missing genotypes for each column
  @inbounds for j in 1:n
    @simd for i in 1:m
      (a1, a2) = A[i, j]
      if isnan(a1, a2)
        missings_by_person[i] += 1
        missings_by_snp[j] += 1
      else
        maf[j] += convert(Float64, a1 + a2) # accumulate A2 allele count
      end
    end
    maf[j] /= 2.0(m - missings_by_snp[j]) # A2 allele frequency
    minor_allele[j] = maf[j] > 0.5
    if minor_allele[j]  # A1 is the minor allele
      maf[j] = 1.0 - maf[j]
    end
  end
  return maf, minor_allele, missings_by_snp, missings_by_person
end

"""
Compute summary statistics of a SnpVector.

# Output:
* `maf` - minor allele frequency
* `minor_allele` - indicate the minor allele is A1 (`true`) or A2 (`false`)
* `missings` - number of missing genotypes
"""
function summarize(A::SnpLike{1})
  m = length(A)
  maf = 0.0                # minor allele frequency
  missings = 0             # no. missing genotypes
  @inbounds @simd for i in 1:m
    (a1, a2) = A[i]
    if isnan(a1, a2)
      missings += 1
    else
      maf += convert(Float64, a1 + a2)
    end
  end
  maf /= 2.0(m - missings) # A2 allele frequency
  minor_allele = maf > 0.5
  if minor_allele          # A1 is the minor allele
    maf = 1.0 - maf
  end
  return maf, minor_allele, missings
end # function summarize

#---------------------------------------------------------------------------#
# Filtering
#---------------------------------------------------------------------------#

"""
    filter(A[, min_success_rate_per_person, min_success_rate_per_person, maxiters])

Filter a SnpMatrix by genotyping success rate.

# Input
- `A`: a SnpMatrix.
- `min_success_rate_per_snp`: threshold for SNP genotyping success rate.
- `min_success_rate_per_person`: threshold for person genotyping success rate.
- `maxiters`: maximum number of filtering iterations.

# Output
- `snp_index`: BitVector indicating remaining SNPs.
- `person_index`: BitVector indicating remaining people.
"""
function filter(
  A::SnpLike{2},
  min_success_rate_per_snp::Float64 = 0.98,
  min_success_rate_per_person::Float64 = 0.98,
  maxiters::Int = 3
  )

  n, p = size(A)
  snp_index = trues(p)
  person_index = trues(n)
  for r in 1:maxiters
    # summary statistics of remaining people/SNPs
    storage = summarize(view(A, person_index, snp_index))
    if (maximum(storage[3]) / countnz(person_index) <
      1.0 - min_success_rate_per_snp) &&
      (maximum(storage[4]) / countnz(snp_index) <
      1.0 - min_success_rate_per_person)
      # all remaining SNPs and people pass success rate threshold
      break
    else
      # some remaining SNPs and people still below success rate threshold
      snp_index[snp_index] .&= (storage[3] / countnz(person_index)
        .< 1.0 - min_success_rate_per_snp)
      person_index[person_index] .&= (storage[4] / countnz(snp_index)
        .< 1.0 - min_success_rate_per_person)
      r == maxiters && "maxiters reached; some SNPs and/or people may still have
        genotyping success rates below threshold."
    end
  end
  return snp_index, person_index
end

#---------------------------------------------------------------------------#
# Empirical kinship matrices
#---------------------------------------------------------------------------#

"""
    grm(A; method=:GRM, memory_limit=2^30, maf_threshold=0.01)

Compute empirical kinship matrix from a SnpMatrix. Missing genotypes are imputed
on the fly according to minor allele frequencies.

# Input
- `A`: a SnpArray

# Optional Arguments
- `method`: `:GRM` (default) or `:MoM`
- `memory_limit`: cap for memory usage in bytes, default `2^30` (2GB)
- `maf_threshold`: SNPs with MAF `<maf_threshold` are excluded, default 0.01
"""
function grm(
  A::SnpLike{2};
  method::Symbol = :GRM,
  memory_limit::Real = 2.0^30,
  maf_threshold::Real = 0.01
  )
  maf, = summarize(A)
  A_filtered = view(A, :, maf .≥ maf_threshold)
  if method == :GRM
    return _grm(A_filtered, memory_limit)
  elseif method == :MoM
    return _mom(A_filtered, memory_limit)
  end
end # function grm

function _grm(A::SnpLike{2}, memory_limit::Real = 2.0^30)
  n, p = size(A)
  Φ = zeros(n, n)
  if 8.0n * p < memory_limit
    snpchunk = convert(Matrix{Float64}, A; model = :additive, impute = true,
      center = true, scale = true)
    BLAS.syrk!('U', 'N', 0.5 / p, snpchunk, 1.0, Φ)
  else
    # chunsize is chosen to have intermediate matrix taking upto memory_limit bytes
    chunksize = ceil(Int, memory_limit / 8.0n)
    snpchunk = zeros(n, chunksize)
    for chunk in 1:floor(Int, p / chunksize)
      J = ((chunk - 1) * chunksize + 1):(chunk * chunksize)
      copy!(snpchunk, view(A, :, J); model = :additive,
        impute = true, center = true, scale = true)
      BLAS.syrk!('U', 'N', 0.5 / p, snpchunk, 1.0, Φ)
    end
    # last chunk
    J = (p - rem(p, chunksize) + 1):p
    if length(J) > 0
      snpchunk = convert(Matrix{Float64}, view(A, :, J); model = :additive,
        impute = true, center = true, scale = true)
      BLAS.syrk!('U', 'N', 0.5 / p, snpchunk, 1.0, Φ)
    end
  end
  # copy to lower triangular part
  LinAlg.copytri!(Φ, 'U')
  return Φ
end # function _grm

function _mom(A::SnpLike{2}, memory_limit::Real = 2.0^30)
  n, p = size(A)
  Φ = zeros(n, n)
  if 8.0n * p < memory_limit # take no more than 1GB memory
    snpchunk = convert(Matrix{Float64}, A; model = :additive, impute = true)
    @inbounds @simd for i in eachindex(snpchunk)
      snpchunk[i] -= 1.0
    end
    BLAS.syrk!('U', 'N', 0.5, snpchunk, 1.0, Φ)
  else
    # chunsize is chosen to have intermediate matrix taking upto memory_limit bytes
    chunksize = ceil(Int, memory_limit / 8.0n)
    snpchunk = zeros(n, chunksize)
    for chunk in 1:floor(Int, p / chunksize)
      J = ((chunk - 1) * chunksize + 1):chunk * chunksize
      copy!(snpchunk, view(A, :, J); model = :additive, impute = true)
      @inbounds @simd for i in eachindex(snpchunk)
        snpchunk[i] -= 1.0
      end
      BLAS.syrk!('U', 'N', 0.5, snpchunk, 1.0, Φ)
    end
    # last chunk
    J = (p - rem(p, chunksize) + 1):p
    if length(J) > 0
      snpchunk = convert(Matrix{Float64}, view(A, :, J);
        model = :additive, impute = true)
      @inbounds @simd for i in eachindex(snpchunk)
        snpchunk[i] -= 1.0
      end
      BLAS.syrk!('U', 'N', 0.5, snpchunk, 1.0, Φ)
    end
  end
  # shift and scale elements of Φ
  c = 0.0
  maf = summarize(A)[1]
  @inbounds for j in eachindex(maf)
    c += maf[j]^2 + (1.0 - maf[j])^2
  end
  a, b = 0.5p - c, p - c
  @inbounds @simd for j in 1:n
    for i in 1:j
      Φ[i, j] += a
      Φ[i, j] /= b
    end
  end
  # copy to lower triangular part
  LinAlg.copytri!(Φ, 'U')
  return Φ
end # function _mom

#---------------------------------------------------------------------------#
# Principal components
#---------------------------------------------------------------------------#

"""
Principal component analysis of SNP data.

# Arguments
* `A`: n-by-p SnpArray.
* `pcs`: number of principal components. Default is 6.

# Output
* `pcscore`: n-by-pcs matrix of principal component scores, or the top `pcs`
  eigen-SNPs.
* `pcloading`: p-by-pcs matrix. Each column is the principal loadings.
* `pcvariance`: princial variances, equivalent to the top `pcs` eigenvalues of
  the sample covariance matrix.

# TODO: maket it work for Integer type matrix
"""
function pca{T <: AbstractFloat}(A::SnpLike{2}, pcs::Integer = 6,
  t::Type{Matrix{T}} = Matrix{Float64})
  n, p = size(A)
  # memory-mapped genotype matrix G, centered and scaled
  G = Mmap.mmap(t, (n, p))
  copy!(G, A; model = :additive, impute = true, center = true, scale = true)
  Gmap = LinearMap(G)
  pcvariance, pcloading, = eigs(Gmap' * Gmap, nev = pcs)
  # make first entry of each eigenvector nonnegative for identifiability
  identify!(pcloading)
  pcscore = G * pcloading
  # scale by n-1 to obtain eigenvalues of the covariance matrix G'G / (n - 1)
  pcvariance ./= n - 1
  return pcscore, pcloading, pcvariance
end # function pca

function pca_sp{T <: Real, TI}(A::SnpLike{2}, pcs::Integer = 6,
  t::Type{SparseMatrixCSC{T, TI}} = SparseMatrixCSC{Float64, Int})
  n, p = size(A)
  # genotype matrix *not* centered or scaled
  G = convert(t, A; model = :additive, impute = true)
  # center and scale
  maf,   = summarize(A)
  center = zeros(T, p)
  weight = zeros(T, p)
  @inbounds @simd for i in 1:p
      center[i] = 2maf[i]
      weight[i] = maf[i] == 0? 1 : 1 / sqrt(2maf[i] * (1 - maf[i]))
  end
  # linear map defined by standardized genotype matrix
  tmpv   = zeros(T, p) # pre-allocate space for intermediate vector
  Gcsv!  = (output, v) ->  Acs_mul_B!(output, G, v, center, weight, tmpv)
  Gcstw! = (output, w) -> Acsc_mul_B!(output, G, w, center, weight)
  Gs     = LinearMap{T}(Gcsv!, Gcstw!, n, p; ismutating=true)
  # PCs
  pcvariance, pcloading, = eigs(Gs' * Gs, nev = pcs)
  # make first entry of each eigenvector nonneagtive for identifiability
  identify!(pcloading)
  # pcscore = Gs * pcloading
  pcscore = zeros(eltype(center), n, pcs)
  Acs_mul_B!(pcscore, G, pcloading, center, weight)
  # scale by n-1 to obtain eigenvalues of the covariance matrix G'G / (n - 1)
  @inbounds @simd for i in 1:pcs
    pcvariance[i] = pcvariance[i] / (n - 1)
  end
  return pcscore, pcloading, pcvariance
end # function pca_sp

"""
Standardized A (centered by `center` and scaled by `weight`) multiply B.
"""
function Acs_mul_B!{T, N}(output::AbstractArray{T, N}, A::AbstractMatrix,
  B::AbstractArray{T, N}, center::Vector{T}, weight::Vector{T},
  storage::AbstractArray{T, N} = similar(B))
  A_mul_B!(storage, Diagonal(weight), B)
  A_mul_B!(output, A, storage)
  @inbounds for j in 1:size(output, 2)
    shift = dot(center, view(storage, :, j))
    @simd for i in 1:size(output, 1)
      output[i, j] -= shift
    end
  end
  output
end # function Acs_mul_B

"""
Transpose of standardized matrix A (centered by `center` and scaled by `weight`)
multiply B.
"""
function Acsc_mul_B!(output::AbstractVector, A::AbstractMatrix,
  b::AbstractVector, center::AbstractVector, weight::AbstractVector)
  Ac_mul_B!(output, A, b)
  shift = sum(b)
  @inbounds @simd for i in eachindex(output)
    output[i] -= shift * center[i]
    output[i] *= weight[i]
  end
  output
end # function Acsc_mul_B!

"""
Make first entry of each column nonnegative. This is for better identifibility
of an orthogonal matrix such as from eigen- or singular value decomposition.
"""
function identify!{T <: AbstractFloat}(A::Matrix{T})
  m, n = size(A)
  zeroT = zero(eltype(A))
  @inbounds for j in 1:n
    if A[1, j] < zeroT
      @simd for i in 1:m
        A[i, j] = - A[i, j]
      end
    end
  end
end # function identify!

#---------------------------------------------------------------------------#
# SnpData type implementation
#---------------------------------------------------------------------------#
include("snpdata.jl")

#---------------------------------------------------------------------------#
# HaplotypeArray type implementation
#---------------------------------------------------------------------------#
include("haplotypearray.jl")

end # module SnpArrays
