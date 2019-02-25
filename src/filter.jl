"""
    SnpArrays.filter(s[, min_success_rate_per_row, min_success_rate_per_col, min_maf, min_hwe_pval, maxiters])

Filter a SnpArray according to genotyping success rate, minor allele frequencies, 
and Hardy-Weinberg test.

# Input
- `s`: a SnpArray.
- `min_success_rate_per_row`: threshold for SNP genotyping success rate. Default 0.98.
- `min_success_rate_per_col`: threshold for person genotyping success rate. Default 0.98.
- `min_maf`: minimum minor allele frequency. Default 0.01.
- `min_hwe_pval`: minimum p-value for Hardy-Weinberg test. Default 0 (not filter HWE).
- `maxiters`: maximum number of filtering iterations. Default is 5.

# Output
- `rmask`: BitVector indicating rows after filtering.
- `cmask`: BitVector indicating cols after filtering.
"""
function filter(
    s::SnpArray, 
    min_success_rate_per_row::Real = 0.98,
    min_success_rate_per_col::Real = 0.98,
    min_maf::Real = 0.01,
    min_hwe_pval::Real = 0,
    maxiters::Integer = 5)
    m, n = size(s)
    # row-wise counts of missing genotypes
    rmissings = zeros(Int, m)
    # columnwise counts: n00, missing, n01, n11
    cc = zeros(Int, 4, n)
    rmask, cmask = trues(m), trues(n)
    rmissrate, cmissrate = 1 - min_success_rate_per_row, 1 - min_success_rate_per_col
    for iter in 1:maxiters
        # accumulate row and column counts
        fill!(rmissings, 0)
        fill!(cc, 0)
        @inbounds for j in 1:n
            cmask[j] || continue
            for i in 1:m
                rmask[i] || continue
                sij = s[i, j]
                cc[sij + 1, j] += 1
                sij == 0x01 && (rmissings[i] += 1)
            end
        end
        rows, cols = count(rmask), count(cmask)
        rmisses, cmisses = rmissrate * rows, cmissrate * cols
        @inbounds for j in 1:n
            cmask[j] || continue
            # if too many missing genotypes, filter out
            if cc[2, j] > rmisses
                cmask[j] = false; continue
            end
            # if maf too low, filter out
            maf = (cc[3, j] + 2cc[4, j]) / 2(cc[1, j] + cc[3, j] + cc[4, j])
            maf = maf < 0.5 ? maf : 1 - maf
            if maf < min_maf
                cmask[j] = false; continue
            end
            # if HWE p-value too low, filter out
            if min_hwe_pval > 0 && hwe(cc[1, j], cc[3, j], cc[4, j]) < min_hwe_pval
                cmask[j] = false; continue
            end
        end
        @inbounds for i in 1:m
            rmask[i] = rmask[i] && rmissings[i] < cmisses
        end
        count(cmask) == cols && count(rmask) == rows && break
        iter == maxiters && @warn("success rate not satisfied; consider increase maxiters")
    end
    rmask, cmask
end

"""
    SnpArrays.filter(src, rowinds, colinds; des = src * ".filtered")

Filter `src` Plink files according to row indices `rowinds` and column indices 
`colinds` and write to a new set of Plink files `des`.

# Input
- `src`: source Plink file name without suffix ".bed", ".fam" or ".bim".
- `rowinds`: row indices.
- `colinds`: column indices.

# Keyword arguments
- `des`: output Plink file name; default is `src * ".filtered"`.
"""
function filter(
    src::AbstractString, 
    rowinds::AbstractVector{<:Integer},
    colinds::AbstractVector{<:Integer};
    des::AbstractString = src * ".filtered")
    # check source plink files
    isfile(src * ".bed") || throw(ArgumentError("$src.bed file not found"))
    isfile(src * ".bim") || throw(ArgumentError("$src.bim file not found"))
    isfile(src * ".fam") || throw(ArgumentError("$src.fam file not found"))
    # create row and column masks
    if eltype(rowinds) == Bool
        rmask = rowinds
    else
        rmask = falses(countlines(src * ".fam"))
        rmask[rowinds] .= true
    end
    if eltype(colinds) == Bool
        cmask = colinds
    else
        cmask = falses(countlines(src * ".bim"))
        cmask[colinds] .= true
    end
    m, n = count(rmask), count(cmask)
    bfsrc = SnpArray(src * ".bed")
    # write filtered bed file
    open(des * ".bed", "w+") do io
        write(io, 0x1b6c)
        write(io, 0x01)
        write(io, Matrix{UInt8}(undef, (m + 3) >> 2, n))
    end
    bfdes = SnpArray(des * ".bed", m, "r+")
    bfdes .= @view bfsrc[rmask, cmask]
    # write filtered fam file
    open(des * ".fam", "w") do io
        for (i, line) in enumerate(eachline(src * ".fam"))
            rmask[i] && println(io, line)
        end
    end
    # write filtered bim file
    open(des * ".bim", "w") do io
        for (j, line) in enumerate(eachline(src * ".bim"))
            cmask[j] && println(io, line)
        end
    end
    # output SnpArray
    bfdes
end

"""
    hwe(n00, n01, n11)

Hardy-Weinberg equilibrium test. `n00`, `n01`, `n11` are counts of homozygotes 
and heterozygoes respectively. Output is the p-value of type Float64.
"""
function hwe(n00::Integer, n01::Integer, n11::Integer)
    n = n00 + n01 + n11
    n == 0 && return 1.0
    p0 = (n01 + 2n00) / 2n
    (p0 ≤ 0.0 || p0 ≥ 1.0) && return 1.0
    p1 = 1 - p0
    # Pearson's Chi-squared test
    e00 = n * p0 * p0
    e01 = 2n * p0 * p1
    e11 = n * p1 * p1
    ts = (n00 - e00)^2 / e00 + (n01 - e01)^2 / e01 + (n11 - e11)^2 / e11
    pval = ccdf(Chi(1), ts)
    # TODO Fisher exact test
    return pval
end