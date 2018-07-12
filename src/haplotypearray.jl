#module HaplotypeArrays

export HaplotypeArray, convert_haplotype

struct HaplotypeArray{N} <: AbstractArray{NTuple{2, Bool}, N}
  A1::BitArray{N}
  A2::BitArray{N}
end

#---------------------------------------------------------------------------#
# Constructors
#---------------------------------------------------------------------------#

# construct HaplotypeArray from SnpArray
HaplotypeArray(s::SnpArray) = HaplotypeArray(s.A1, s.A2)

# construct SnpArray from HaplotypeArray
function SnpArray(h::HaplotypeArray)
  s = SnpArray(h.A1, h.A2)
  # NaN (true, false) should be (false, true)
  @simd for i in eachindex(s)
    if isnan(s[i]); s[i] = (false, true); end
  end
  return s
end

function HaplotypeArray(dims...)
  HaplotypeArray(falses(dims), falses(dims))
end #

# HaplotypeArray or a view of a HaplotypeArray
@compat HaplotypeLike = Union{HaplotypeArray{N}, SubArray{NTuple{2, Bool}, N, H}} where {N, H<:HaplotypeArray}
@compat HaplotypeMatrix = HaplotypeArray{2}
@compat HaplotypeVector = HaplotypeArray{1}

#---------------------------------------------------------------------------# methods
# Julia docs on methods required for AbstractArray:
# http://docs.julialang.org/en/release-0.4/manual/interfaces/#man-interfaces-abstractarray

Base.size(A::HaplotypeArray)                 = size(A.A1)
Base.size(A::HaplotypeArray, d::Int)         = size(A.A1, d)
Base.ndims(A::HaplotypeArray)                = ndims(A.A1)
Base.endof(A::HaplotypeArray)                = length(A)
Base.eltype(A::HaplotypeArray)               = NTuple{2, Bool}
@compat Base.IndexStyle(::Type{HaplotypeArray})      = Base.IndexLinear()

@inline function Base.getindex(A::HaplotypeArray, i::Int)
  (getindex(A.A1, i), getindex(A.A2, i))
end

@inline function Base.getindex(A::HaplotypeArray, i::Int, j::Int)
  (getindex(A.A1, i, j), getindex(A.A2, i, j))
end

@inline function Base.setindex!(A::HaplotypeArray, v::NTuple{2, Bool}, i::Int)
  setindex!(A.A1, v[1], i), setindex!(A.A2, v[2], i)
end

@inline function Base.setindex!(A::HaplotypeArray, v::NTuple{2, Bool}, i::Int, j::Int)
  setindex!(A.A1, v[1], i, j), setindex!(A.A2, v[2], i, j)
end

function Base.similar(A::HaplotypeArray, dims::Dims)
  HaplotypeArray(BitArray(dims), BitArray(dims))
end # function Base.similar

function Base.similar(A::HaplotypeArray, ::Type{NTuple{2, Bool}}, dims::Dims)
  similar(A, dims)
end # function Base.similar

#---------------------------------------------------------------------------#
# Code for missing genotypes
#---------------------------------------------------------------------------#

Base.isnan(A::HaplotypeLike{N}) where {N} = falses(size(A)) 

#---------------------------------------------------------------------------#
# Convert and copy
#---------------------------------------------------------------------------#

"""
Create a HaplotypeArray with all A1 alleles.
"""
function Base.convert(t::Type{HaplotypeArray}, dims...)
  HaplotypeArray(falses(dims), falses(dims))
end # function Base.convert

"""
Convert a two-bit genotype to a real number (minor allele count) according to
specified SNP model, ignoring the missing genotype code. `minor_allele=true`
indicates `A1` is the minor allele; `minor_allele=false` indicates `A2` is the
minor allele.
"""
@inline function convert_haplotype(t::Type{T}, a::NTuple{2, Bool},
  minor_allele::Bool, model::Symbol = :additive) where {T <: Real}
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
  return b
end # function convert_haplotype

"""
Convert a haplotype matrix to a numeric matrix of minor allele counts according to
specified SNP model.
"""
function Base.convert(t::Type{Array{T, N}}, A::HaplotypeLike{N};
  model::Symbol = :additive, center::Bool = false, scale::Bool = false) where {T <: Real, N}
  B = similar(A, T)
  copy!(B, A; model = model, center = center, scale = scale)
end # function Base.convert

function Base.copy!(
  B::Array{T, N},
  A::Union{HaplotypeLike{1}, HaplotypeLike{2}};
  model::Symbol = :additive,
  impute::Bool = false,
  center::Bool = false,
  scale::Bool = false
  ) where {T <: Real, N}

  @assert size(B) == size(A) "Dimensions do not match"
  m, n = size(A, 1), size(A, 2)
  # convert column by column
  @inbounds for j in 1:n
    # first pass: find minor allele and its frequency
    if N == 1
      maf, minor_allele = summarize(A)
    elseif N == 2
      maf, minor_allele = summarize(view(A, :, j))
    end
    # second pass: convert, center, scale
    ct = convert(T, 2.0maf)
    wt = convert(T, maf == 0.0 ? 1.0 : 1.0 / √(2.0maf * (1.0 - maf)))
    @simd for i in 1:m
      (a1, a2) = A[i, j]
      B[i, j] = convert_haplotype(T, (a1, a2), minor_allele, model)
      if center; B[i, j] -= ct; end
      if scale; B[i, j] *= wt; end
    end
  end
  return B
end # function Base.copy!

"""
Convert a SNP matrix to a numeric matrix of minor allele counts according to
specified SNP model.
"""
function Base.convert(
  t::Type{SparseMatrixCSC{T, TI}},
  A::HaplotypeLike{2};
  model::Symbol = :additive
  ) where {T <: Real, TI <: Integer}

  m, n = size(A)
  # prepare sparese matrix data structure
  rowval = TI[]
  nzval = T[]
  colptr = zeros(TI, n + 1)
  # convert column by column
  zeroT = zero(T)
  v = zeroT
  @inbounds for j in 1:n
    colptr[j] = convert(TI, length(nzval) + 1)
    # first pass: find minor allele and its frequency
    maf, minor_allele = summarize(view(A, :, j))
    # second pass: impute, convert
    for i in 1:m
      (a1, a2) = A[i, j]
      v = convert_haplotype(T, (a1, a2), minor_allele, model)
      if v ≠ zeroT
        push!(rowval, convert(TI, i)), push!(nzval, v)
      end
    end # i
  end # j
  colptr[n + 1] = length(nzval) + 1
  return SparseMatrixCSC(m, n, colptr, rowval, nzval)
end # function Base.convert

#---------------------------------------------------------------------------#
# Summary statistics
#---------------------------------------------------------------------------#

"""
Compute summary statistics of a HaplotypeMatrix.

# Output:
* `maf` - minor allele frequency of each SNP
* `minor_allele` - indicate the minor allele is A1 (`true`) or A2 (`false`)
"""
function summarize(A::HaplotypeLike{2})

  m, n = size(A)
  maf = zeros(Float64, n)             # minor allele frequencies for each column
  minor_allele = trues(n)             # true->A1 is the minor allele
  @inbounds for j in 1:n
    @simd for i in 1:m
      (a1, a2) = A[i, j]
      maf[j] += convert(Float64, a1 + a2) # accumulate A2 allele count
    end
    maf[j] /= 2m # A2 allele frequency
    minor_allele[j] = maf[j] > 0.5
    if minor_allele[j]  # A1 is the minor allele
      maf[j] = 1.0 - maf[j]
    end
  end
  return maf, minor_allele
end

"""
Compute summary statistics of a HaplotypeVector.

# Output:
* `maf` - minor allele frequency
* `minor_allele` - indicate the minor allele is A1 (`true`) or A2 (`false`)
"""
function summarize(A::HaplotypeLike{1})

  m = length(A)
  maf = 0.0                # minor allele frequency
  missings = 0             # no. missing genotypes
  @inbounds @simd for i in 1:m
    (a1, a2) = A[i]
    maf += convert(Float64, a1 + a2)
  end
  maf /= 2m # A2 allele frequency
  minor_allele = maf > 0.5
  if minor_allele          # A1 is the minor allele
    maf = 1.0 - maf
  end
  return maf, minor_allele
end # function summarize

#end # HaplotypeArrays module
