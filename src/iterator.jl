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
    maf_vector = similar(s.snparray, Float64)
    maf!(maf_vector, s.snparray)
    return maf_vector[snpindex.index]
end

function hwepval(s::SnpData, snpindex::SnpArrayIndex)
    genotypes = s.snparray[:,snpindex.index]

    n00 = sum(genotypes .== 0) 
    n01 = sum(genotypes .== 1) 
    n11 = sum(genotypes .== 2) # / length(genotypes)

   pval = hwe(n00,n01,n11)
   return pval
    
    # function hwe()
    # filter.jl line 188

end

function alt_dosages!(arr::AbstractArray{T}, s::SnpData, snpindex::SnpArrayIndex; use_genotype::Bool=false) where T <: Real
    allele1, allele2 = alleles(s, snpindex)
    
    for i in 1:size(arr, 1)
        for j in 1:size(arr, 2)
            if s.snparray[i, snpindex.index] == 0x01
                arr[i, j] = 0.5  # Heterozygous dosage
            elseif s.snparray[i, snpindex.index] == 0x02
                arr[i, j] = 1.0  # Homozygous dosage for allele2
            else
                arr[i, j] = 0.0  # Homozygous dosage for allele1 or missing
            end
            
            if use_genotype
                if s.snparray[i, snpindex.index] == 0x01
                    arr[i, j] = 1  # Heterozygous
                elseif s.snparray[i, snpindex.index] == 0x02
                    arr[i, j] = 2  # Homozygous for allele2
                else
                    arr[i, j] = missing  # Missing genotype or homozygous for allele1
                end
            end
        end
    end
    
    return arr
end



function alt_genotypes!(arr::AbstractArray{T}, s::SnpData, snpindex::SnpArrayIndex) where T <: Integer
    allele1, allele2 = alleles(s, snpindex)
    
    @inbounds for i in eachindex(arr)
        if s.snparray[i, snpindex.index] == 0x01
            arr[i] = 1  # Heterozygous
        elseif s.snparray[i, snpindex.index] == 0x02
            arr[i] = 2  # Homozygous for allele2
        else
            arr[i] = missing  # Missing genotype or homozygous for allele1
        end
    end
    return arr
end
