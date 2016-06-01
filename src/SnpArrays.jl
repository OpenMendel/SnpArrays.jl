module SnpArrays

import IterativeSolvers: MatrixFcn
export estimatesize, grm, pca, pca_sp, randgeno, SnpArray, SnpData, summarize

type SnpArray{N} <: AbstractArray{NTuple{2, Bool}, N}
  A1::BitArray{N}
  A2::BitArray{N}
end

# SnpArray or a view of a SnpArray
typealias SnpLike{N} Union{SnpArray{N}, SubArray{NTuple{2, Bool}, N}}
typealias SnpMatrix SnpArray{2}
typealias SnpVector SnpArray{1}

"""
Construct a SnpArray from an array of A1 allele counts {0, 1, 2}.
"""
function SnpArray(mac::AbstractArray)
    T = eltype(mac)
    SnpArray(mac .> one(T), mac .> zero(T))
end

"""
Construct a SnpArray from Plink binary files (`.bed`, `.bim`, `.fam`).
"""
function SnpArray(plinkFile::AbstractString)
  plinkBedfile = string(plinkFile, ".bed")
  plinkBimfile = string(plinkFile, ".bim")
  plinkFamfile = string(plinkFile, ".fam")
  # dimensions
  nper = countlines(plinkFamfile)
  nsnp = countlines(plinkBimfile)
  # read binary genotype data from bed file
  fid = open(plinkBedfile, "r")
  # decide BED file version and snp/individual-major
  bedheader = read(fid, UInt8, 3)
  # PLINK coding (genotype->bits): A1/A1->00, A1/A2->01, A2/A2->11, missing->10
  # http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml
  A1 = BitArray(nper, nsnp)
  A2 = BitArray(nper, nsnp)
  if bits(bedheader[1]) == "01101100" && bits(bedheader[2]) == "00011011"
    # v1.0 BED file
    if bits(bedheader[3]) == "00000001"
      # SNP-major
      plinkbits = Mmap.mmap(fid, BitArray, (2, 4ceil(Int, 0.25nper), nsnp), 3)
      A1 = copy!(A1, slice(plinkbits, 1, 1:nper, :))
      A2 = copy!(A2, slice(plinkbits, 2, 1:nper, :))
    else
      # individual-major
      snpbits = Mmap.mmap(fid, BitArray, (2, 4ceil(Int, 0.25nsnp), nper), 3)
      A1 = copy!(A1, slice(plinkbits, 1, 1:nper, :)')
      A2 = copy!(A2, slice(plinkbits, 2, 1:nper, :)')
    end
  else
    # TODO: v0.99 BED file: individual-major
    error("openmendel:plinkformat\n",
          "v0.99 BED file found!",
          "Transform to v1.0 BED file using PLINK")
  end
  close(fid)
  SnpArray(A1, A2)
end

#---------------------------------------------------------------------------# methods
# Julia docs on methods required for AbstractArray:
# http://docs.julialang.org/en/release-0.4/manual/interfaces/#man-interfaces-abstractarray
Base.size(A::SnpArray)                 = size(A.A1)
Base.size(A::SnpArray, d::Int)         = size(A.A1, d)
Base.ndims(A::SnpArray)                = ndims(A.A1)
Base.endof(A::SnpArray)                = length(A)
Base.eltype(A::SnpArray)               = NTuple{2, Bool}
Base.linearindexing(::Type{SnpArray})  = Base.LinearFast()

function Base.getindex(A::SnpArray, i::Int)
  (getindex(A.A1, i), getindex(A.A2, i))
end

function Base.getindex(A::SnpArray, i::Int, j::Int)
  (getindex(A.A1, i, j), getindex(A.A2, i, j))
end

function Base.setindex!(A::SnpArray, v::NTuple{2, Bool}, i::Int)
  setindex!(A.A1, v[1], i), setindex!(A.A2, v[2], i)
end

function Base.setindex!(A::SnpArray, v::NTuple{2, Bool}, i::Int, j::Int)
  setindex!(A.A1, v[1], i, j), setindex!(A.A2, v[2], i, j)
end

function Base.setindex!(A::SnpArray, v::Real, i::Int)
  if isnan(v)
    setindex!(A, (true, false), i)
  else
    setindex!(A, (v > one(eltype(v)), v > zero(eltype(v))), i)
  end
end

function Base.setindex!(A::SnpArray, v::Real, i::Int, j::Int)
  if isnan(v)
    setindex!(A, (true, false), i, j)
  else
    setindex!(A, (v > one(eltype(v)), v > zero(eltype(v))), i, j)
  end
end

function Base.similar(A::SnpArray, T::Real, dims::Dims)
  zeros(T, dims)
end

function Base.similar(A::SnpArray)
  SnpArray(similar(A.A1), similar(A.A2))
end

function Base.similar(A::SnpArray, ::NTuple{2, Bool})
  similar(A)
end

function Base.similar(A::SnpArray, dims::Dims)
  SnpArray(BitArray(dims), BitArray(dims))
end

function Base.similar(A::SnpArray, ::NTuple{2, Bool}, dims::Dims)
  SnpArray(BitArray(dims), BitArray(dims))
end

# missing code is 10 = (true, false)
Base.isnan(a1::Bool, a2::Bool) = a1 & !a2
Base.isnan(a::Tuple{Bool, Bool}) = Base.isnan(a[1], a[2])
function Base.isnan{N}(A::SnpLike{N})
  b = BitArray(size(A))
  @inbounds @simd for i in eachindex(A)
    b[i] = Base.isnan(A[i])
  end
  return b
end

"""
Create a SnpArray with all A1 alleles.
"""
function Base.convert(t::Type{SnpArray}, dims...)
  SnpArray(falses(dims), falses(dims))
end

"""
Convert a two-bit genotype to a real number according to specified SNP model.
Missing genotype is converted to `NaN`. `minor_allele=true` indicates `A1` is
the minor allele; `minor_allele=false` indicates `A2` is the minor allele.
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
end

"""
Convert a SNP matrix to a numeric matrix according to specified SNP model. If
`impute == true`, missing entries are imputed according to (column) minor allele
frequencies.
"""
function Base.convert{T <: Real, N}(t::Type{Array{T, N}}, A::SnpLike{N};
  model::Symbol = :additive, impute::Bool = false, center::Bool = false,
  scale::Bool = false)
  B = similar(A, T)
  copy!(B, A; model = model, impute = impute, center = center, scale = scale)
end

function Base.copy!{T <: Real, N}(B::Array{T, N}, A::SnpLike{N};
  model::Symbol = :additive, impute::Bool = false, center::Bool = false,
  scale::Bool = false)
  @assert size(B) == size(A) "Dimensions do not match"
  if ndims(A) == 1
    m, n = length(A), 1
  elseif ndims(A) == 2
    m, n = size(A)
  end
  # convert column by column
  @inbounds for j in 1:n
    # first pass: find minor allele and its frequency
    maf, minor_allele, = summarize(sub(A, :, j))
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
      if center
        B[i, j] -= ct
      end
      if scale
        B[i, j] *= wt
      end
    end
  end
  return B
end

"""
Convert a SNP matrix to a sparse matrix according to specified SNP model.
If `impute=false`, missing genotypes are ignored (translated to 0).
If `impute=true`, missing genotypes are imputed on the fly according
to the minor allele frequencies.
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
    maf, minor_allele, = summarize(sub(A, :, j))
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
end

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
end

"""
Generate a genotype according to minor allele frequency. `minor_allele` indicates
the minor allele is A1 (`true`) or A2 (`false`).
"""
function randgeno{T <: AbstractFloat}(maf::T, minor_allele::Bool)
  minor_allele ? randgeno(maf) : randgeno(1.0 - maf)
end

"""
Generate a SnpVector according to minor allele frequency. `minor_allele` indicates
the minor allele is A1 (`true`) or A2 (`false`).
"""
function randgeno{T <: AbstractFloat}(n::Int, maf::T, minor_allele::Bool)
  s = SnpArray(n)
  @inbounds @simd for i in 1:n
    s[i] = randgeno(maf, minor_allele)
  end
  return s
end

"""
Generate a SnpMatrix according to minor allele frequencies `maf`. `minor_allele`
vector indicates the minor alleles are A1 (`true`) or A2 (`false`).
"""
function randgeno{T <: AbstractFloat}(m::Int, n::Int, maf::Vector{T},
  minor_allele::BitVector)
  @assert length(maf) == n "length of maf should be n"
  @assert length(minor_allele) == n "length of minor_allele should be n"
  s = SnpArray(m, n)
  @inbounds @simd for j in 1:n
    for i in 1:m
      s[i, j] = randgeno(maf[j], minor_allele[j])
    end
  end
  return s
end

function randgeno{T <: AbstractFloat}(n::Tuple{Int,Int}, maf::Vector{T},
  minor_allele::BitArray{1})
  randgeno(n..., maf, minor_allele)
end


"""
Compute summary statistics of a SnpMatrix.

# Output:
* `maf` - minor allele frequency of each SNP
* `minor_allele` - indicate the minor allele is A1 (`true`) or A2 (`false`)
*  `missings_by_person` - number of missing genotypes for each person
*  `missings_by_snp` - number of missing genotypes for each SNP
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
end

"""
Compute empirical kinship matrix from a SnpMatrix. Missing genotypes are imputed
on the fly according to minor allele frequencies.
"""
function grm(A::SnpLike{2}; method::Symbol = :GRM)
  if method == :GRM
    return _grm(A::SnpLike{2})
  elseif method == :MoM
    return _mom(A::SnpLike{2})
  end
end

function _grm(A::SnpLike{2})
  n, p = size(A)
  Φ = zeros(n, n)
  memory_limit = 2.0^30 # 1GB memory usage limit
  if 8.0n * p < memory_limit
    snpchunk = convert(Matrix{Float64}, A; model = :additive, impute = true,
      center = true, scale = true)
    BLAS.syrk!('U', 'N', 0.5 / p, snpchunk, 1.0, Φ)
  else
    # chunsize is chosen to have intermediate matrix taking upto 1GB memory
    chunksize = ceil(Int, memory_limit / 8.0n)
    snpchunk = zeros(n, chunksize)
    for chunk in 1:floor(Int, p / chunksize)
      J = ((chunk - 1) * chunksize + 1):(chunk * chunksize)
      copy!(snpchunk, sub(A, :, J); model = :additive,
        impute = true, center = true, scale = true)
      BLAS.syrk!('U', 'N', 0.5 / p, snpchunk, 1.0, Φ)
    end
    # last chunk
    J = (p - rem(p, chunksize) + 1):p
    if length(J) > 0
      snpchunk = convert(Matrix{Float64}, sub(A, :, J); model = :additive,
        impute = true, center = true, scale = true)
      BLAS.syrk!('U', 'N', 0.5 / p, snpchunk, 1.0, Φ)
    end
  end
  # copy to lower triangular part
  LinAlg.copytri!(Φ, 'U')
  return Φ
end

function _mom(A::SnpLike{2})
  n, p = size(A)
  memory_limit = 2.0^30 # 1GB memory usage limit
  Φ = zeros(n, n)
  if 8.0n * p < memory_limit # take no more than 1GB memory
    snpchunk = convert(Matrix{Float64}, A; model = :additive, impute = true)
    @inbounds @simd for i in eachindex(snpchunk)
      snpchunk[i] -= 1.0
    end
    BLAS.syrk!('U', 'N', 0.5, snpchunk, 1.0, Φ)
  else
    # chunsize is chosen to have intermediate matrix taking upto 1GB memory
    chunksize = ceil(Int, memory_limit / 8.0n)
    snpchunk = zeros(n, chunksize)
    for chunk in 1:floor(Int, p / chunksize)
      J = ((chunk - 1) * chunksize + 1):chunk * chunksize
      copy!(snpchunk, sub(A, :, J); model = :additive, impute = true)
      @inbounds @simd for i in eachindex(snpchunk)
        snpchunk[i] -= 1.0
      end
      BLAS.syrk!('U', 'N', 0.5, snpchunk, 1.0, Φ)
    end
    # last chunk
    J = (p - rem(p, chunksize) + 1):p
    if length(J) > 0
      snpchunk = convert(Matrix{Float64}, sub(A, :, J);
        model = :additive, impute = true)
      @inbounds @simd for i in eachindex(snpchunk)
        snpchunk[i] -= 1.0
      end
      BLAS.syrk!('U', 'N', 0.5, snpchunk, 1.0, Φ)
    end
  end
  # shift and scale elements of Φ
  c = 0.0
  maf, = summarize(A)
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
end

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
function pca{T <: AbstractFloat}(A::SnpLike{2}, pcs::Int = 6,
  t::Type{Matrix{T}} = Matrix{Float64})
  n, p = size(A)
  # memory-mapped genotype matrix G, centered and scaled
  G = Mmap.mmap(t, (n, p))
  copy!(G, A; model = :additive, impute = true, center = true, scale = true)
  # partial SVD
  _, pcvariance, pcloading = svds(G, nsv = pcs)
  # make first entry of each eigenvector nonneagtive for identifiability
  identify!(pcloading)
  pcscore = G * pcloading
  # square singular values and scale by n - 1 to get eigenvalues of the
  # covariance matrix
  @inbounds @simd for i in 1:pcs
    pcvariance[i] = pcvariance[i] * pcvariance[i] / (n - 1)
  end
  return pcscore, pcloading, pcvariance
end

function pca_sp{T <: Real, TI}(A::SnpLike{2}, pcs::Int = 6,
  t::Type{SparseMatrixCSC{T, TI}} = SparseMatrixCSC{Float64, Int})
  n, p = size(A)
  # genotype matrix *not* centered or scaled
  G = convert(t, A; model = :additive, impute = true)
  # center and scale
  maf, = summarize(A)
  center = 2.0maf
  weight = map((x) -> x == 0.0 ? 1.0 : 1.0 / √(2.0x * (1.0 - x)), maf)
  # standardized genotype matrix
  tmpv = zeros(eltype(center), n) # pre-allocate space for intermediate vector
  Gs = MatrixFcn{eltype(center)}(p, p,
    (output, v) -> AcstAcs_mul_B!(output, G, v, center, weight, tmpv))
  # PCs
  pcvariance, pcloading = eigs(Gs, nev = pcs)
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
end

"""
Gram matrix Acs'Acs multiply B, where Acs is the centered and scaled A.
"""
function AcstAcs_mul_B!{T}(output::Vector{T}, A::AbstractMatrix,
  b::Vector{T}, center::Vector{T}, weight::Vector{T},
  tmp::Vector{T} = zeros(eltype(b), size(A, 1)))
  # output is used as a temporary vector here
  @inbounds @simd for i in eachindex(output)
    output[i] = weight[i] * b[i]
  end
  A_mul_B!(tmp, A, output)
  shift = dot(center, output)
  @inbounds @simd for i in eachindex(tmp)
    tmp[i] -= shift
  end
  # now output is the real output
  Ac_mul_B!(output, A, tmp)
  @inbounds @simd for i in eachindex(output)
    output[i] *= weight[i]
  end
  output
end

"""
Cheating: to bypass the isssym() error thrown by arpack.jl
"""
Base.issym{T}(fcn::MatrixFcn{T}) = true

"""
Standardized A (centered by `center` and scaled by `weight`) multiply B.
"""
function Acs_mul_B!{T, N}(output::AbstractArray{T, N}, A::AbstractMatrix,
  B::AbstractArray{T, N}, center::Vector{T}, weight::Vector{T},
  tmp::AbstractArray{T, N} = zeros(eltype(B), size(B)))
  @inbounds @simd for j in 1:size(B, 2)
    for i in 1:size(B, 1)
      tmp[i, j] = weight[i] * B[i, j]
    end
  end
  A_mul_B!(output, A, tmp)
  @inbounds for j in 1:size(output, 2)
    shift = dot(center, sub(tmp, :, j))
    @simd for i in 1:size(output, 1)
      output[i, j] -= shift
    end
  end
  output
end

"""
Transpose of standardized matrix A (centered by `center` and scaled by `weight`)
multiply B.
"""
function Acsc_mul_B!(output::AbstractVector, A::AbstractMatrix,
  b::AbstractVector, center::AbstractVector, weight::AbstractVector,
  tmpvec::AbstractVector = zeros(eltype(A), size(A, 2)))
  Ac_mul_B!(tmpvec, A, b)
  shift = sum(b)
  @inbounds @simd for i in eachindex(tmpvec)
    tmpvec[i] = (tmpvec[i] - shift * center[i]) * weight[i]
  end
  output
end

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
end

"""
Estimate the memory usage for storing SNP data with `n` individuals,
`p` SNPs, and average minor allele frequency `maf`.
"""
function estimatesize{T <: AbstractFloat}(n::Int, p::Int,
  maf::T, t::Type)
  if t <: DenseMatrix
    storage = convert(Float64, sizeof(eltype(t)) * n * p)
  elseif t <: AbstractSparseMatrix
    storage = convert(Float64,
      (sizeof(t.parameters[1]) + sizeof(t.parameters[2])) *
      n * p * maf * (2.0 - maf) + sizeof(t.parameters[2]) * (p + 1))
  end
  storage
end

"""
Type for SNP and person information.
"""
type SnpData
  people::Int                       # number of rows (individuals)
  snps::Int                         # number of columns (snps)
  personid::Vector{AbstractString}  # names of individuals
  snpchr::Vector{AbstractString}    # snp chromosome
  snppos::Vector{AbstractString}    # snp position
  snpid::Vector{AbstractString}     # ids of snps
  maf::Vector{Float64}              # minor allele frequencies
  minor_allele::BitVector            # bit array designating the minor allele
  snpmatrix::SnpLike{2}             # matrix of genotypes or haplotypes
  missings_per_person::Vector{Int}  # number of missing genotypes per person
  missings_per_snp::Vector{Int}     # number of missing genotypes per snp
end # end type

"""
Construct a SnpData type from a PLINK file.
"""
function SnpData(plink_file::AbstractString)

  # load snp info
  plink_bim_file = string(plink_file, ".bim")
  snp_info = readdlm(plink_bim_file, AbstractString)
  snpchr = snp_info[:, 1]
  snpid = snp_info[:, 2]
  snppos = snp_info[:, 4]
  #allele1 = snp_info[:, 5]
  #allele2 = snp_info[:, 6]

  # load person info
  plink_fam_file = string(plink_file, ".fam")
  person_info = readdlm(plink_fam_file, AbstractString)
  personid = person_info[:, 2]

  # load snp array matrix
  snpmatrix = SnpArray(plink_file)
  maf, minor_allele, missings_per_snp, missings_per_person = summarize(snpmatrix)
  people, snps = size(snpmatrix)

  # construct SnpData unfiltered
  snp_data = SnpData(people, snps, personid, snpchr, snppos, snpid,
    maf, minor_allele, snpmatrix, missings_per_person, missings_per_snp)

end

end # module
