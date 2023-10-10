mutable struct SnpArrayIterator <: VariantIterator
   snpdata::SnpData
end

mutable struct SnpArrayIndex <: Variant
    index::Int
end 

@inline function Base.eltype(::Type{<:VariantIterator})
    SnpArrayIndex
end

function Base.iterate(itr::SnpArrayIterator, state=1)
    if state <= 0
        throw(BoundsError(itr, state))
    end
    if state > size(itr.snpdata.snp_info,2)
        return nothing
    else
        index = SnpArrayIndex(state)
        # don't need view just return the snparray index and the next state 
        # tuple of length 2 
        state = state + 1
        return (index, state)
    end
end

@inline function Base.length(itr::SnpArrayIterator)
    return size(itr.snpdata.snp_info, 2)
end

function iterator(s::SnpData; snpindex::SnpArrayIndex)
    iterator = SnpArrayIterator(s)
end

function chrom(s::SnpData, snpindex::SnpArrayIndex)::String
    result = s.snp_info[snpindex.index,:chromosome]
    return result
end

function pos(s::SnpData, snpindex::SnpArrayIndex)::Int
    result = s.snp_info[snpindex.index,:position]
    return result
end


function rsid(s::SnpData, snpindex::SnpArrayIndex)::String
    result = s.snp_info[snpindex.index,:snpid]
    return result
end

#SnpData subtype of Genetic Data

function alleles(s::SnpData, snpindex::SnpArrayIndex)::Vector{String}
    allele1 = s.snp_info[snpindex.index,:allele1]
    allele2 = s.snp_info[snpindex.index,:allele2]
    return [allele1, allele2]
end

function alt_allele(s::SnpData, snpindex::SnpArrayIndex)::String
    alt = s.snp_info[snpindex.index,:allele2]
    return alt
end

function ref_allele(s::SnpData, snpindex::SnpArrayIndex)::String
    ref = s.snp_info[snpindex.index,:allele1]
    return ref
end

function maf(s::SnpData, snpindex::SnpArrayIndex)
    array = s.snparray
    return maf(array)
end

function hwepval(s::SnpData, snpindex::SnpArrayIndex)'
# Extract the genotype data for a specific SNP using snpindex
    genotypes = s.snparray[snpindex, :]

    # Calculate the observed genotype frequencies
    observed_hom_ref = sum(genotypes .== 0) / length(genotypes)
    observed_het = sum(genotypes .== 1) / length(genotypes)
    observed_hom_alt = sum(genotypes .== 2) / length(genotypes)

    # Calculate the expected genotype frequencies under HWE
    total_alleles = 2 * s.people
    expected_hom_ref = (observed_hom_ref + observed_het / 2)
    expected_het = observed_het
    expected_hom_alt = (observed_hom_alt + observed_het / 2)
    
end