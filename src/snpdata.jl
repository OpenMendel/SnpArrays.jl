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
    src::AbstractString
end

function SnpData(plink_file::AbstractString)
    
    # load snp info
    plink_bim_file = string(plink_file, ".bim")
    snp_info = categorical!(
                   CSV.read(plink_bim_file, delim='\t', header=SNP_INFO_KEYS, 
                        types=[String, String, Float64, Int, String, String]),
                   [:allele1, :allele2])
    
    # load person info
    plink_fam_file = string(plink_file, ".fam")
    person_info = convert(DataFrame, readdlm(plink_fam_file, AbstractString))
    rename!(person_info, f => t for (f, t) = zip(names(person_info),
                PERSON_INFO_KEYS))

    # load snp array 
    snparray = SnpArray(string(plink_file, ".bed"))
    people, snps = size(snparray)
    
    SnpData(people, snps, snparray, snp_info, person_info, plink_file)
end

"""
    show(io::IO, x::SnpData)

Simple string representation of SnpData
"""
function Base.show(io::IO, x::SnpData)
    print(io, "SnpData(people: $(x.people), snps: $(x.snps),\n" *
        join([
                "snp_info: \n$((join(split(string(first(x.snp_info, 6)), "\n")[2:end], "\n")))",
                x.snps > 6 ? "...," : ",",
                "person_info: \n$((join(split(string(first(x.person_info, 6)), "\n")[2:end], "\n")))",
                x.people > 6 ? "...," : ",",
                "src: $(x.src)"], 
                "\n") *
        "\n)")
end


"""
    filter(s::SnpData, rowinds::AbstractVector{<:Integer}, colinds::AbstractVector{<:Integer}; des::AbstractString)
    filter(s; des::AbstractString, f_person::Function, f_snp::Function)

Filter `s` according to `f_person` and `f_snp`. The resultiing plink files are saved at `des.[bed|bim|fam]`.
"""
function filter(s::SnpData, rowinds::AbstractVector{<:Integer}, colinds::AbstractVector{<:Integer}; des::AbstractString = s.src * ".filtered")
    snps = sum(colinds)
    people = sum(rowinds)
    snparray = SnpArrays.filter(s.src, rowinds, colinds; des=des)
    snp_info = s.snp_info[colinds, :]
    person_info = s.person_info[rowinds, :]
    SnpData(people, snps, snparray, snp_info, person_info, des)
end
function filter(s::SnpData; des::AbstractString = s.src * ".filtered", f_person::Function = (x -> true), f_snp::Function = (x -> true))
    colinds = collect(f_snp(r)::Bool for r in eachrow(s.snp_info)) 
    rowinds = collect(f_person(r)::Bool for r in eachrow(s.person_info))
    SnpArrays.filter(s::SnpData, rowinds::Vector{Bool}, colinds::Vector{Bool}; des = des::AbstractString)
end
filter(s::AbstractString; des::AbstractString = s.src * ".filtered", 
    f_person::Function = (x -> true), f_snp::Function = (x -> true)) = SnpArrays.filter(SnpData(s); des = des, f_person = f_person, f_snp = f_snp)


"""
    split_plink(s, key; prefix)

Split data `s` according to chromosome, sex, or phenotype. Returns a dictionary of splitted data.
"""
function split_plink(s::SnpData, key::Symbol = :chromosome; prefix = s.src * string(key))
    key in [:chromosome, :sex, :phenotype] || throw(ArgumentError("key should be one of :chromosome, :sex, or :phenotype"))
    r = Dict{AbstractString, SnpData}()
    if key == :chromosome
        for chr in unique(s.snp_info[key])
            r[chr] = SnpArrays.filter(s; des = prefix * chr, f_snp = x -> (x[key] == chr))
        end
        r
    else
        for val in unique(s.person_info[key])
            r[val] = SnpArrays.filter(s; des = prefix * string(val), f_person = x -> (x[key] == val))
        end
        r
    end
end
split_plink(src::AbstractString, key::Symbol = :chromosome; prefix = s.src * string(key)) = split_plink(SnpData(src), key; prefix = prefix)


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
    SnpData(people, snps, snparray, snp_info, person_info, "")
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
    # isfile(bimfile) && error("($bimfile) already exists.")
    # isfile(bedfile) && error("($bedfile) already exists.")
    # isfile(famfile) && error("($famfile) already exists.")
    
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
    SnpData(snpdata.people, snpdata.snps, snpdata.snparray, snpdata.snp_info, snpdata.person_info, filename)
end

"""
    compresse_plink(plink_file, format[="gz"], outfile[=plink_file])

Compress a set of Plink files to `gz` (default), `zlib` or `zz` format.
"""
function compress_plink(
    plink_file::AbstractString, 
    format::AbstractString="gz", 
    outfile::AbstractString=plink_file
    )
    if format == "gz"
        stream = GzipCompressorStream
    elseif format == "zlib"
        stream = ZlibCompressorStream
    elseif format == "zz"
        stream = DeflateCompressorStream
    else 
        throw(ArgumentError("compressed format should be gz or zlib or zz"))
    end
    for suffix in [".bed", ".fam", ".bim"]
        fname = outfile * suffix
        if isfile(fname)
            open(fname) do input
                open(stream, fname * "." * format, "w") do output
                    write(output, input)
                end
            end
        end
    end
end

"""
    decompresse_plink(plink_compressed, outfile[=plink_compressed])

Decompress a set of compressed Plink files to plain Plink format.
"""
function decompress_plink(
    plink_compressed::AbstractString, 
    format::AbstractString="gz",
    outfile::AbstractString=plink_compressed
    )
    if format == "gz"
        stream = GzipDecompressorStream
    elseif format == "zlib"
        stream = ZlibDecompressorStream
    elseif format == "zz"
        stream = DeflateDecompressorStream
    else 
        throw(ArgumentError("compressed format should be gz or zlib or zz"))
    end
    for suffix in [".bed", ".fam", ".bim"]
        iname = plink_compressed * suffix * "." * format
        oname = outfile * suffix
        if isfile(iname)
            open(oname, "w") do output
                open(stream, iname) do input
                    write(output, input)
                end
            end
        end
    end
end
