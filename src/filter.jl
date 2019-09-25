"""
    SnpArrays.filter(s)

Filter a SnpArray according to genotyping success rate, minor allele frequencies, 
and/or Hardy-Weinberg test.

# Input
- `s`: a SnpArray or Plink file name without the bim, fam, bed suffix.

# Keyword argument
- `min_success_rate_per_row`: Threshold for SNP genotyping success rate. Default 0.98. 
- `min_success_rate_per_col`: Threshold for person genotyping success rate. Default 0.98. 
- `min_maf`: Minimum minor allele frequency. Default 0.01.
- `min_hwe_pval`: Minimum p-value for Hardy-Weinberg test. Default 0 (not filter HWE).
- `maxiters`: Maximum number of filtering iterations. Default is 5.

# Output
- `rmask`: BitVector indicating rows after filtering.
- `cmask`: BitVector indicating columns after filtering.
"""
function filter(
    s::SnpArray;
    min_success_rate_per_row::Real = 0.98,
    min_success_rate_per_col::Real = 0.98,
    min_maf::Real = 0.01,
    min_hwe_pval::Real = 0,
    maxiters::Integer = 5)
    m, n = size(s)
    # row-wise counts of missing genotypes
    rmissings = zeros(Int, m)
    # create temporary missings to add to rmissings if chromosome is counted (not filtered)
    temprmissings = zeros(Int, m)
    # columnwise counts: n00, missing, n01, n11
    cc = zeros(Int, 4)
    rmask, cmask = trues(m), trues(n)
    rmissrate, cmissrate = 1 - min_success_rate_per_row, 1 - min_success_rate_per_col
    for iter in 1:maxiters
        # number of remaining rows
        rows = count(rmask)
        # maximum allowed missing genotypes each row and col
        cmisses = cmissrate * rows
        fill!(rmissings, 0)
        @inbounds for j in 1:n
            cmask[j] || continue
            # accumulate row and column counts
            fill!(cc, 0)
            fill!(temprmissings, 0)
            for i in 1:m
                rmask[i] || continue
                sij = s[i, j]
                cc[sij + 1] += 1
                sij == 0x01 && (temprmissings[i] += 1)
            end
            # if too many missing genotypes, filter out
            if cc[2] > cmisses
                cmask[j] = false; continue
            end
            # if maf too low, filter out
            maf = (cc[3] + 2cc[4]) / 2(cc[1] + cc[3] + cc[4])
            maf = maf ≤ 0.5 ? maf : 1 - maf
            if maf < min_maf
                cmask[j] = false; continue
            end
            # if HWE p-value too low, filter out
            if min_hwe_pval > 0 && hwe(cc[1], cc[3], cc[4]) < min_hwe_pval
                cmask[j] = false; continue
            end
            # if snp passes all filters, then count the row missings
            rmissings += temprmissings
        end
        # number of remaining cols 
        cols = count(cmask)
        # maximum allowed missing genotypes each row 
        rmisses = rmissrate * cols
        # filter rows/samples
        @inbounds for i in 1:m
            rmask[i] = rmask[i] && rmissings[i] < rmisses
        end
        # if no change in filter results, done
        count(cmask) == cols && count(rmask) == rows && break
        iter == maxiters && @warn("success rate not satisfied; consider increase maxiters")
    end
    rmask, cmask
end

# function filter(src::AbstractString; des::AbstractString = src * ".filtered", kwargs...)
#     dirname = splitdir(src)[1]
#     srcbedfile = readdir(glob"src.bed", dirname)[1]
#     srcfamfile = readdir(glob"src.fam", dirname)[1]
#     srcm = makestream(srcfamfile) do stream
#         countlines(stream)
#     end
#     s = SnpArray(srcbedfile, srcm)
#     rowmask, colmask = SnpArrays.filter(s; kwargs...)
#     SnpArrays.filter(src, rowmask, colmask; des=des)
# end

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
    # source bed, fam, bim file names
    dirname, filename = splitdir(src)
    srcbedfile = glob(filename * ".bed", dirname)[1]
    srcbimfile = glob(filename * ".bim", dirname)[1]
    srcfamfile = glob(filename * ".fam", dirname)[1]
    # check source plink files
    isfile(srcbedfile) || throw(ArgumentError("$src.bed file not found"))
    isfile(srcbimfile) || throw(ArgumentError("$src.bim file not found"))
    isfile(srcfamfile) || throw(ArgumentError("$src.fam file not found"))
    # destination bed, fam, bim file names
    desbedfile = replace(srcbedfile, src => des)
    desbimfile = replace(srcbimfile, src => des)
    desfamfile = replace(srcfamfile, src => des)
    # numbers of samples and SNPs in src
    srcm = makestream(srcfamfile) do stream
        countlines(stream)
    end
    srcn = makestream(srcbimfile) do stream
        countlines(stream)
    end
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
    desm, desn = count(rmask), count(cmask)
    # write filtered bed file
    bfsrc = SnpArray(srcbedfile)
    makestream(desbedfile, "w+") do io
        write(io, 0x1b6c)
        write(io, 0x01)
        write(io, Matrix{UInt8}(undef, (desm + 3) >> 2, desn))
    end
    bfdes  = SnpArray(desbedfile, desm, "r+")
    bfdes .= @view bfsrc[rmask, cmask]
    # write filtered fam file
    makestream(desfamfile, "w") do io
        for (i, line) in enumerate(eachline(srcfamfile))
            rmask[i] && println(io, line)
        end
    end
    # write filtered bim file
    makestream(desbimfile, "w") do io
        for (j, line) in enumerate(eachline(srcbimfile))
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
    pval = ccdf(Chisq(1), ts)
    # TODO Fisher exact test
    return pval
end


"""
    indexin_sorted(v, w)

Same as `indexin_general(v, w)` but assumes `v` and `w` are sorted.
"""
function indexin_sorted(v::Union{AbstractArray, Tuple}, w)
    viter, witer = keys(v), eachindex(w)
    vind, wmask = eltype(viter)[], falses(length(w))
    vy, wy = iterate(viter), iterate(witer)
    if vy === nothing || wy === nothing
        return vind
    end
    viteri, i = vy
    witerj, j = wy
    @inbounds begin
        vi, wj = v[viteri], w[witerj]
        while true
            if isless(vi, wj)
                # advance vi
                vy = iterate(viter, i)
                if vy === nothing
                    break
                end
                viteri, i = vy
                vi        = v[viteri]
            elseif isless(wj, vi)
                # advance wj
                wy = iterate(witer, j)
                if wy === nothing
                    break
                end
                witerj, j = wy
                wj        = w[witerj]
            else
                push!(vind, viteri)
                wmask[witerj] = true
                # advance vi
                vy = iterate(viter, i)
                if vy === nothing
                    break
                end
                viteri, i = vy
                vi        = v[viteri]
                # advance wj
                wy = iterate(witer, j)
                if wy === nothing
                    break
                end
                witerj, j = wy
                wj        = w[witerj]
            end
        end
    end
    return vind, wmask
end


"""
    indexin_general(v, w)

Returns an index vector `vind` and `wmask` such that `v[vind]` is the subset 
of `v` that appear in `w` and `wmask` is a bitvector such that `v[vind] .== w[wmask]`. 
Repeated matches in `v` will be ignored. Only the first match is kept.

# Examples
```jldoctest
julia> a = ['b', 'c', 'd'];

julia> b = ['a', 'b', 'c'];

julia> SnpArrays.indexin_general(a, b)
([1, 2], Bool[false, true, true])
```
"""
function indexin_general(v::Union{AbstractArray, Tuple}, w)
    if issorted(v) && issorted(w)
        return indexin_sorted(v, w)
    else
        Iv, Iw = sortperm(v), sortperm(w)
        vind, wmask = indexin_sorted(v[Iv], w[Iw])
        invpermute!(wmask, Iw)
        if v[Iv[vind]] == w[wmask]
            return Iv[vind], wmask
        else
            return sort!(Iv[vind]), wmask
        end
    end
end

"""
    SnpArrays.filter(srcbedfile, srcbimfile, srcfamfile, rowinds, colinds; des = src * ".filtered")

Filter Plink files  with .gz format or differently named bim and bed files according to row indices 
`rowinds` and column indices `colinds` and write to a new set of Plink files `des`.

# Input
- `srcbedfile`: bed file name with suffix such as .bed or .bed.gz.
- `srcbimfile`: bed file name with suffix such as .bim or .bim.gz.
- `srcfamfile`: bed file name with suffix such as .fam or .fam.gz.
- `rowinds`: row indices.
- `colinds`: column indices.

# Keyword arguments
- `des`: output Plink file name; default is `src * ".filtered"`.
"""
function filter(
    srcbedfile::AbstractString,
    srcbimfile::AbstractString,
    srcfamfile::AbstractString, 
    rowinds::AbstractVector{<:Integer},
    colinds::AbstractVector{<:Integer};
    des::Union{Nothing, AbstractString} = nothing)
    srcbed = split(srcbedfile, ".bed")[1]
    srcbim = split(srcbimfile, ".bim")[1]
    srcfam = split(srcfamfile, ".fam")[1]
    # check source plink files
    isfile(srcbedfile) || throw(ArgumentError("$srcbedfile file not found"))
    isfile(srcbimfile) || throw(ArgumentError("$srcbimfile file not found"))
    isfile(srcfamfile) || throw(ArgumentError("$srcfamfile file not found"))
    # destination bed, fam, bim file names
    if des == nothing
        desbedfile = srcbed * ".filtered.bed"
        desbimfile = srcbim * ".filtered.bim"
        desfamfile = srcfam * ".filtered.fam"
    else
        desbedfile = des * ".bed"
        desbimfile = des * ".bim"
        desfamfile = des * ".fam"
    end
    # numbers of samples and SNPs in src
    srcm = makestream(srcfamfile) do stream
        countlines(stream)
    end
    srcn = makestream(srcbimfile) do stream
        countlines(stream)
    end
    # create row and column masks
    if eltype(rowinds) == Bool
        rmask = rowinds
    else
        rmask = falses(srcm)
        rmask[rowinds] .= true
    end
    if eltype(colinds) == Bool
        cmask = colinds
    else
        cmask = falses(srcn)
        cmask[colinds] .= true
    end
    desm, desn = count(rmask), count(cmask)
    # write filtered bed file
    bfsrc = SnpArray(srcbedfile, srcm)
    makestream(desbedfile, "w+") do io
        write(io, 0x1b6c)
        write(io, 0x01)
        write(io, Matrix{UInt8}(undef, (desm + 3) >> 2, desn))
    end
    bfdes  = SnpArray(desbedfile, desm, "r+")
    bfdes .= @view bfsrc[rmask, cmask]
    # write filtered fam file
    makestream(desfamfile, "w") do io
        if endswith(srcfamfile, ".gz")
            gzio = GzipDecompressorStream(open(srcfamfile, "r"))
            for (i, line) in enumerate(eachline(gzio))
                rmask[i] && println(io, line)
            end
        else
            for (i, line) in enumerate(eachline(srcfamfile))
                rmask[i] && println(io, line)
            end
        end
    end
    # write filtered bim file
    makestream(desbimfile, "w") do io
        if endswith(srcbimfile, ".gz")
            gzio = GzipDecompressorStream(open(srcbimfile, "r"))
            for (j, line) in enumerate(eachline(gzio))
                cmask[j] && println(io, line)
            end
        else
            for (j, line) in enumerate(eachline(srcbimfile))
                cmask[j] && println(io, line)
            end
        end
    end
    # output SnpArray
    bfdes
end