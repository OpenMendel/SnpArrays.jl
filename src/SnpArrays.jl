module SnpArrays

export SnpArray, SnpArrayBM, summarysnps, grm, impute!, summarize

type SnpArray{N} <: AbstractArray{NTuple{2, Bool}, N}
  A1::BitArray{N}
  A2::BitArray{N}
end

"""
For bench marking purpose
"""
type SnpArrayBM{N} <: AbstractArray{NTuple{2, Bool}, N}
  A::BitArray{N}
end

"""
Construct a SnpArray from an array of minor allele counts {0, 1, 2}.
# TODO: make it robust with NaN entries
"""
function SnpArray(mac::AbstractArray)
    T = eltype(mac)
    SnpArray(mac .> zero(T), mac .> one(T))
end


"""
Construct a SnpArray from Plink binary files.
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
  # PLINK coding (bits->genotype): 00->aa, 01->Aa, 10->Missing, 11->AA
  # where a is the minor allele. We code it as bitwise inverse of Plink.
  # Our coding: 11->aa, 10->Aa, 01->Missing, 00->AA.
  if bits(bedheader[1]) == "01101100" && bits(bedheader[2]) == "00011011"
    # v1.0 BED file
    if bits(bedheader[3]) == "00000001"
      # SNP-major
      snpbits = Mmap.mmap(fid, BitArray{3}, (2, 4ceil(Int, 0.25nper), nsnp), 3; grow=false)
      A1 = !slice(snpbits, 1, 1:nper, :)
      A2 = !slice(snpbits, 2, 1:nper, :)
    else
      # individual-major
      snpbits = Mmap.mmap(fid, BitArray{3}, 
        (2, 4ceil(Int, 0.25nsnp), nper), 3)
      A1 = !slice(snpbits, 1, 1:nsnp, :)'
      A2 = !slice(snpbits, 2, 1:nsnp, :)'
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

"""
Construct a SnpArrayBM from Plink binary files.
"""
function SnpArrayBM(plinkFile::AbstractString)
  plinkBedfile = string(plinkFile, ".bed")
  plinkBimfile = string(plinkFile, ".bim")
  plinkFamfile = string(plinkFile, ".fam")
  nper = countlines(plinkFamfile)
  nsnp = countlines(plinkBimfile)
  fid = open(plinkBedfile, "r")
  bedheader = read(fid, UInt8, 3)
  if bits(bedheader[1]) == "01101100" && bits(bedheader[2]) == "00011011"
    if bits(bedheader[3]) == "00000001"
      snpbits = Mmap.mmap(fid, BitArray{3}, (2, 4ceil(Int, 0.25nper), nsnp),
        3; grow=false, shared=false)
      A = !slice(snpbits, :, 1:nper, :)
    else
      snpbits = Mmap.mmap(fid, BitArray{3}, (2, 4ceil(Int, 0.25nsnp), nper),
        3; grow=false, shared=false)
      A = !slice(snpbits, :, 1:nsnp, :)'
    end
  else
    error("openmendel:plinkformat\n",
          "v0.99 BED file found!",
          "Transform to v1.0 BED file using PLINK")
  end
  close(fid)
  SnpArrayBM(A)
end

# SnpArray or a view of a SnpArray
typealias SnpLike{N} Union{ SnpArray{N}, SnpArrayBM{N}, 
    SubArray{NTuple{2, Bool}, N}}
typealias SnpMatrix SnpArray{2}
typealias SnpVector SnpArray{1}

#---------------------------------------------------------------------------# methods
# Julia docs on methods required for AbstractArray:
# http://docs.julialang.org/en/release-0.4/manual/interfaces/#man-interfaces-abstractarray
Base.size(A::SnpArray)                 = size(A.A1)
Base.size(A::SnpArray, d::Int)         = size(A.A1, d)

"""
For SnpArrayBM
"""
Base.size(A::SnpArrayBM)                 = size(A.A)

"""
For SnpArrayBM
"""
Base.size(A::SnpArrayBM, d::Int)         = size(A.A, d)

Base.ndims(A::SnpArrayBM)                = ndims(A.A)
Base.ndims(A::SnpArray)                = ndims(A.A1)
Base.endof(A::SnpArray)                = length(A)
Base.eltype(A::SnpArray)               = NTuple{2, Bool}
Base.eltype(A::SnpArrayBM)               = NTuple{2, Bool}
Base.linearindexing(::Type{SnpArray})  = Base.LinearSlow()

"""
For SnpArrayBM
"""
Base.linearindexing(::Type{SnpArrayBM})  = Base.LinearSlow()

function Base.getindex(A::SnpArray, i::Int)
  (getindex(A.A1, i), getindex(A.A2, i))
end

function Base.getindex(A::SnpArray, i::Int, j::Int)
  (getindex(A.A1, i, j), getindex(A.A2, i, j))
end

"""
Get genotype for individual i and SNP j for SnpArrayBM
"""
function Base.getindex(SA::SnpArrayBM, i::Int, j::Int)
  #(getindex(SA.A, 1, i, j), getindex(SA.A, 2, i, j))
  (SA.A[1,i,j], SA.A[2,i,j])
end

function Base.setindex!(A::SnpArray, v::NTuple{2, Bool}, i::Int)
  setindex!(A.A1, v[1], i), setindex!(A.A2, v[2], i)
end

function Base.setindex!(A::SnpArray, v::NTuple{2, Bool}, i::Int, j::Int)
  setindex!(A.A1, v[1], i, j), setindex!(A.A2, v[2], i, j)
end

"""
Set genotype for individual i and SNP j for SnpArrayBM
"""
function Base.setindex!(SA::SnpArrayBM, v::NTuple{2, Bool}, i::Int, j::Int)
  setindex!(SA.A, v[1], 1, i, j), setindex!(SA.A, v[2], 2, i, j)
end

function Base.setindex!(A::SnpArray, v::Real, i::Int)
  if isnan(v)
    setindex!(A, (false, true), i)
  else
    setindex!(A, (v > zero(eltype(v)), v > one(eltype(v))), i)
  end
end

function Base.setindex!(A::SnpArray, v::Real, i::Int, j::Int)
  if isnan(v)
    setindex!(A, (false, true), i, j)
  else
    setindex!(A, (v > zero(eltype(v)), v > one(eltype(v))), i, j)
  end
end

function Base.similar(A::SnpArray, T::Real, dims::Dims)
  zeros(T, dims)
end

function Base.similar(A::SnpArrayBM{3}, T::Type)
  _,m,n = size(A)
  zeros(T, (m,n))
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
Convert a two-bit genotype to a real number. Missing genotype is converted to
NaN.
"""
function Base.convert{T<:Real}(t::Type{T}, a::NTuple{2, Bool};
  model::Symbol = :additive)
  b = zero(t)
  if !a[1] & a[2]
    b = convert(t, NaN)
  else
    if model == :additive
      b = convert(T, a[1] + a[2])
    elseif model == :dominant
      b = convert(T, 2(a[1] & a[2]))
    elseif model == :recessive
      b = convert(T, 2a[1])
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

"""
Convert function for SnpArrayBM
"""
function Base.convert{T <: Real, N}(t::Type{Array{T, N}}, A::SnpLike{3};
  model::Symbol = :additive, impute::Bool = false, center::Bool = false,
  scale::Bool = false)
  B = similar(A, T)
  copy!(B, A; model = model, impute = impute, center = center, scale = scale)
end

function Base.copy!{T <: Real, N}(B::Array{T, N}, A::SnpLike;
  model::Symbol = :additive, impute::Bool = false, center::Bool = false,
  scale::Bool = false)
  nanT = convert(T, NaN)
  if ndims(A) == 1
    m, n = length(A), 1
  elseif ndims(A) == 2
    m, n = size(A)
  elseif ndims(A) == 3
    _, m, n = size(A)
  end
  
  @inbounds for j = 1:n
    # first pass: convert and count missing genotypes
    nmisscol = 0 # no. missing entries in column j
    nmialcol = 0  # no. minor alleles in column j
    @simd for i = 1:m
      (a1, a2) = A[i, j]
      if !a1 & a2 # missing value (false, true)
        B[i, j] = nanT
        nmisscol += 1
      else
        nmialcol += a1 + a2
        if model == :additive
          B[i, j] = a1 + a2
        elseif model == :dominant
          B[i, j] = 2(a1 & a2)
        elseif model == :recessive
          B[i, j] = 2a1
        end
      end
    end

    # second pass: impute, center, scale
    if (impute && nmisscol > 0) || center || scale
      maf = nmialcol / 2(m - nmisscol)
      ct = 2.0maf
      wt = ifelse(maf == 0.0, 1.0, 1.0 / sqrt(2.0maf * (1.0 - maf)))
      @inbounds @simd for i = 1:m
        if isnan(B[i, j])
          a1 = rand() < maf
          a2 = rand() < maf
          if model == :additive
            B[i, j] = a1 + a2
          elseif model == :dominant
            B[i, j] = 2(a1 & a2)
          elseif model == :recessive
            B[i, j] = 2a1
          end
        end
        if center
          B[i, j] -= ct
        end
        if scale
          B[i, j] *= wt
        end
      end
    end
  end
  return B
end

"""
Convert a SNP matrix to a sparse matrix according to specified SNP model.
Missing genotypes are ignored.
#TODO: implement imputation.
"""
function Base.convert{T <: Real, TI <: Integer, N}(t::Type{SparseMatrixCSC{T, TI}},
  A::SnpLike{N}; model::Symbol = :additive)
  if ndims(A) == 2
    m, n = size(A)
  elseif ndims(A) == 3
    _, m, n = size(A)
  end
  I, J, V = TI[], TI[], T[]
  zeroT = zero(T)
  v = zero(T)
  @inbounds for j = 1:n, i = 1:m
    (a1, a2) = A[i, j]
    if !a1 & a2 # missing
      continue
    end
    if model == :additive
      v = convert(T, a1 + a2)
    elseif model == :dominant
      v = convert(T, 2(a1 & a2))
    elseif model == :recessive
      v = convert(T, 2a1)
    end
    if v ≠ zeroT
      push!(I, i), push!(J, j), push!(V, v)
    end
  end
  return sparse(I, J, V, m, n)
end


Base.isnan(a::Tuple{Bool, Bool}) = !a[1] & a[2]
Base.isnan(A::SnpArray) = !A.A1 & A.A2
# TODO: make isnan work for view
# function Base.isnan{N}(A::SubArray{Tuple{Bool, Bool}, N, SnpArray{N},
#   Tuple{UnitRange{Int64}, N}, 0})
#   pA = parent(A)
#   pidx = parentindexes(A)
#   !pA.A1[pidx] & pA.A2[pidx]
# end

# TODO:
function impute!(A::SnpLike{2})

end

"""
Computes the number of minor alleles, number of missing genotypes, and minor
allele frequencies (MAF) along each row and column. Calculation of MAF takes
into account of missingness.
"""
function summarysnps(A::SnpLike{2})
  m, n = size(A)
  nmialcol = zeros(Int, n)      # no. minor alleles for each column
  nmisscol = zeros(Int, n)      # no. missing genotypes for each column
  mafcol   = zeros(Float64, n)  # minor allele frequencies for each column
  nmialrow = zeros(Int, m)      # no. minor alleles for each row
  nmissrow = zeros(Int, m)      # no. missing genotypes for each row
  mafrow   = zeros(Float64, m)  # minor allele frequencies for each row
  @inbounds for j = 1:n
    @simd for i = 1:m
      (a1, a2) = A[i, j]
      ism = !a1 & a2
      nmialcol[j] += ifelse(ism, 0, a1 + a2)
      nmisscol[j] += ism
      nmialrow[i] += ifelse(ism, 0, a1 + a2)
      nmissrow[i] += ism
    end
    mafcol[j] = nmialcol[j] / 2(m - nmisscol[j])
  end
  @inbounds @simd for i = 1:m
    mafrow[i] = nmialrow[i] / 2(m - nmissrow[i])
  end
  return nmialcol, nmisscol, mafcol, nmialrow, nmissrow, mafrow
end


"""
Compute summary statistics of a SnpMatrix.
# Output:
* `maf` - minor allele frequency of each SNP
* `minor_allele` - indicate the minor allele is A1 (`true`) or A2 (`false`)
*  `missings_by_person` - number of missing genotypes for each person
*  `missings_by_snp` - number of missing genotypes for each SNP
"""
function summarize(A::SnpLike)
  if ndims(A) == 2
    m, n = size(A)
  elseif ndims(A) == 3
    _,m,n = size(A)
  end
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
Computes empirical kinship matrix from a SnpMatrix.
# TODO: implement MoM method
"""
function grm(A::SnpLike; method::Symbol = :GRM)
  if method == :GRM
    return _grm(A)
  elseif method == :MoM
    return _mom(A)
  end
end

function _grm(A::SnpLike)
  if ndims(A) == 2
    n, p = size(A)
  elseif ndims(A) == 3
    _, n, p = size(A)
  end
  Φ = zeros(n, n)
  memory_limit = 2.0^30 # 1GB memory usage limit
  if 8.0n * p < memory_limit
    snpchunk = convert(Matrix{Float64}, A; model = :additive, impute = true,
      center = true, scale = true)
    BLAS.syrk!('U', 'N', 0.5 / p, snpchunk, 1.0, Φ)
  else
    # chunsize is chosen to have intermediate matrix taking upto 1GB memory
    chunksize = ceil(Int, memory_limit / 8.0n)
    snpchunk = zeros(Float64, (n, chunksize))
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

function _mom(A::SnpLike)
  if ndims(A) == 2
    n, p = size(A)
  elseif ndims(A) == 3
    _, n, p = size(A)
  end
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

end # module


S = SnpArrays
