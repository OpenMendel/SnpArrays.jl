const SNP_INFO_KEYS = [:chromosome, :snpid, :genetic_distance, :position, :allele1, :allele2]
const PERSON_INFO_KEYS = [:fid, :iid, :father, :mother, :sex, :phenotype]

"""
    SnpData
    SnpData(plknm)

Type to store SNP and person information along with the SnpArray.
"""
struct SnpData
    people::Int
    snps::Int
    snparray::SnpArray
    snp_info::DataFrame  
    person_info::DataFrame 
end

function SnpData(plink_file::AbstractString)
    
    # load snp info
    plink_bim_file = string(plink_file, ".bim")
    snp_info = convert(DataFrame, readdlm(plink_bim_file, AbstractString))
    rename!(snp_info, f => t for (f, t) = zip(names(snp_info), 
                SNP_INFO_KEYS))
    snp_info[:genetic_distance] = map(x -> parse(Float64, x), 
        snp_info[:genetic_distance])
    snp_info[:position] = map(x -> parse(Int, x), snp_info[:position])
    snp_info[:allele1] = map(x -> x[1], snp_info[:allele1])
    snp_info[:allele2] = map(x -> x[1], snp_info[:allele2])
    
    # load person info
    plink_fam_file = string(plink_file, ".fam")
    person_info = convert(DataFrame, readdlm(plink_fam_file, AbstractString))
    rename!(person_info, f => t for (f, t) = zip(names(person_info),
                PERSON_INFO_KEYS))

    # load snp array 
    snparray = SnpArray(string(plink_file, ".bed"))
    people, snps = size(snparray)
    
    SnpData(people, snps, snparray, snp_info, person_info)
end

"""
    split_plink(src; prefix)

Split `src` Plink files or SnpData according to chromosome. returns: a dictionary of splitted data keyed by name of chromosome 
"""
function split_plink(src::AbstractString; prefix::AbstractString = src * ".chr.")
    data = SnpData(src)
    r = Dict{AbstractString, SnpData}()
    for chr in unique(data.snp_info[:chromosome])
        ind = (chr .== data.snp_info[:chromosome])
        subarray = filter(src, trues(data.people), ind; des = prefix * chr)
        r[chr] = SnpData(prefix * chr)
    end
    r
end

"""
    merge_plink(d)
    merge_plink(des, d)

Merge the SnpData of the splitted plink files. 
Returns: merged SnpData. If `des` is given, it is written on that destination.
"""
function merge_plink(d::Dict{AbstractString, SnpData})
    ks = sort(collect(keys(d)))
   
    # get person_info
    person_info = d[ks[1]].person_info
    
    # vcat snp_info
    snp_info = vcat([d[k].snp_info for k in ks]...)
    
    # hcat snparray
    data = hcat([d[k].snparray.data for k in ks]...)
    rowcounts = d[ks[1]].snparray.rowcounts
    columncounts = hcat([d[k].snparray.columncounts for k in ks]...)
    snparray = SnpArray(data, rowcounts, columncounts, size(data)[1])

    people, snps = size(person_info,1), size(snp_info, 1)
    SnpData(people, snps, snparray, snp_info, person_info)
end

merge_plink(des::AbstractString, d::Dict{AbstractString, SnpData}) = write_plink(des, merge_plink(d))


"""
    merge_plink(prefix; des::AbstractString = prefix * ".merged")

merge the plink files beginning with `prefix`. 
"""
function merge_plink(prefix::AbstractString; des::AbstractString = prefix * ".merged")
    l = glob(prefix * "*.bed")
    matching_srcs = map(x -> splitext(x)[1], l)
    d = Dict{AbstractString, SnpData}()
    for fn in matching_srcs
        chrsnpdata = SnpData(fn)
        chr = chrsnpdata.snp_info[:chromosome][1]
        @assert all(chr .== chrsnpdata.snp_info[:chromosome]) "Not all chrs are the same in $fn.bim."
        d[chr] = chrsnpdata
    end
    merge_plink(des, d)
end

"""
    write_plink(filename, snpdata)

Write SnpData to Plink bed, bim, fam files. 
Returns: `snpdata`
"""
function write_plink(filename::AbstractString, snpdata::SnpData)
    bimfile = filename * ".bim"
    bedfile = filename * ".bed"
    famfile = filename * ".fam"
    isfile(bimfile) && error("($bimfile) already exists.")
    isfile(bedfile) && error("($bedfile) already exists.")
    isfile(famfile) && error("($famfile) already exists.")
    
    # write bim file
    snp_info = snpdata.snp_info
    writedlm(bimfile, hcat([snp_info[k] for k in SNP_INFO_KEYS]...))
    
    # write fam file
    person_info = snpdata.person_info
    writedlm(famfile, hcat([person_info[k] 
                                for k in PERSON_INFO_KEYS]...))
    
    # write bed file
    open(bedfile, "w") do io
        write(io, 0x1b6c)
        write(io, 0x01)
        write(io, snpdata.snparray.data)
    end
    snpdata
end