const SNP_INFO_KEYS = [:chromosome, :snpid, :genetic_distance, :position, :allele1, :allele2]
const PERSON_INFO_KEYS = [:fid, :iid, :father, :mother, :sex, :phenotype]

"""
    SnpData
    SnpData(plknm; famnm::AbstractString=plinknm*".fam", bimnm::AbstractString=plinknm*".bim" )

Type to store SNP and person information along with the SnpArray.
"""
struct SnpData
    people::Int
    snps::Int
    snparray::SnpArray
    snp_info::DataFrame  
    person_info::DataFrame 
    srcbed::AbstractString
    srcbim::AbstractString
    srcfam::AbstractString
end

function SnpData(plink_file::AbstractString, args...; famnm::AbstractString=plink_file*".fam", bimnm::AbstractString=plink_file*".bim", kwargs...)
    
    # load snp info
    snp_info = open(bimnm) do io
        categorical!(CSV.read(io,  delim='\t', header=SNP_INFO_KEYS, 
        types=[String, String, Float64, Int, String, String]), [:allele1, :allele2])
    end
    
    # load person info
    person_info = convert(DataFrame, readdlm(famnm, AbstractString))
    rename!(person_info, collect(f => t for (f, t) = zip(names(person_info),
                    PERSON_INFO_KEYS)))

    # load snp array 
    snparray = SnpArray(string(plink_file, ".bed"), args...; famnm=famnm, kwargs...)
    people, snps = size(snparray)
    
    SnpData(people, snps, snparray, snp_info, person_info, plink_file*".bed", bimnm, famnm)
end

Base.size(x::SnpData) = Base.size(x.snparray)

"""
    show(io::IO, x::SnpData)

Simple string representation of SnpData
"""
function Base.show(io::IO, x::SnpData)
    print(io, "SnpData(people: $(x.people), snps: $(x.snps),\n" *
        join([
                "snp_info: \n$((join(split(string(first(x.snp_info, 6)), "\n")[2:end], "\n")))",
                x.snps > 6 ? "…," : ",",
                "person_info: \n$((join(split(string(first(x.person_info, 6)), "\n")[2:end], "\n")))",
                x.people > 6 ? "…," : ",",
                "srcbed: $(x.srcbed)\nsrcbim: $(x.srcbim)\nsrcfam: $(x.srcfam)"], 
                "\n") *
        "\n)")
end



"""
    filter(s::SnpData, rowinds::AbstractVector{<:Integer}, colinds::AbstractVector{<:Integer}; des::AbstractString)
    filter(s; des::AbstractString, f_person::Function, f_snp::Function)

Filter `s` according to `f_person` and `f_snp`. The resultiing plink files are saved at `des.[bed|bim|fam]`.
"""
@inline _trueftn(x) = true
function filter(s::SnpData, rowinds::AbstractVector{<:Integer}, colinds::AbstractVector{<:Integer}; des::AbstractString=split(s.srcbed, ".bed")[1] * ".filtered")
    if eltype(rowinds) == Bool
        rmask = rowinds
    else
        rmask = falses(size(s.snparray, 1))
        rmask[rowinds] .= true
    end
    if eltype(colinds) == Bool
        cmask = colinds
    else
        cmask = falses(size(s.snparray, 2))
        cmask[colinds] .= true
    end

    snps = count(cmask)
    people = count(rmask)
    snparray = SnpArrays.filter(s.srcbed, s.srcbim, s.srcfam, rmask, cmask; des=des)
    snp_info = s.snp_info[cmask, :]
    person_info = s.person_info[rmask, :]
    SnpData(people, snps, snparray, snp_info, person_info, 
            des * ".bed", des * ".bim", des * ".fam")
end

function filter(s::SnpData; des::AbstractString = split(s.srcbed, ".bed")[1] * ".filtered", f_person::Function = _trueftn, f_snp::Function = _trueftn)
    f_person == _trueftn && f_snp == _trueftn && @warn "No nontrivial function provided. Just copying."
    colinds = collect(f_snp(r)::Bool for r in eachrow(s.snp_info)) 
    rowinds = collect(f_person(r)::Bool for r in eachrow(s.person_info))
    SnpArrays.filter(s::SnpData, rowinds::Vector{Bool}, colinds::Vector{Bool}; des = des::AbstractString)
end

filter(s::AbstractString; des::AbstractString = split(s.srcbed, ".bed")[1] * ".filtered", 
    f_person::Function = _trueftn, f_snp::Function = _trueftn) = SnpArrays.filter(SnpData(s); des = des, f_person = f_person, f_snp = f_snp)


"""
    split_plink(s, key; prefix)

Split data `s` according to chromosome, sex, or phenotype. Returns a dictionary of splitted data.
"""
function split_plink(s::SnpData, key::Symbol = :chromosome; prefix = split(s.srcbed, ".bed")[1] * string(key))
    key in [:chromosome, :sex, :phenotype] || throw(ArgumentError("key should be one of :chromosome, :sex, or :phenotype"))
    r = Dict{AbstractString, SnpData}()
    if key == :chromosome
        for chr in unique(s.snp_info[!, key])
            r[chr] = SnpArrays.filter(s; des = prefix * chr, f_snp = x -> (x[key] == chr))
        end
        r
    else
        for val in unique(s.person_info[!, key])
            r[val] = SnpArrays.filter(s; des = prefix * string(val), f_person = x -> (x[key] == val))
        end
        r
    end
end
split_plink(src::AbstractString, key::Symbol = :chromosome; prefix = src * string(key)) = split_plink(SnpData(src), key; prefix = prefix)

"""
    vcat(A...;des="tmp_vcat_" * string(vcat_counter))

Concatenate SnpData along dimension 1.
"""
vcat_counter = 1
function Base.vcat(A::SnpData...; des="tmp_vcat_" * string(vcat_counter))
    global vcat_counter
    if des == "tmp_vcat_" * string(vcat_counter)
        vcat_counter = vcat_counter + 1
    end
    
    # vcat person_info
    person_info = vcat([x.person_info for x in A]...)
    
    # get snp_info
    @assert allequal([x.snp_info for x in A]) "snp_info are different"
    snp_info = A[1].snp_info
    
    # vcat snparray
    snparray = vcat([x.snparray for x in A]...; des=des)

    people, snps = size(person_info,1), size(snp_info, 1)

    bimfile = des * ".bim"
    famfile = des * ".fam"

    writedlm(bimfile, hcat([snp_info[!, k] for k in SNP_INFO_KEYS]...))
    writedlm(famfile, hcat([person_info[!, k] 
                                for k in PERSON_INFO_KEYS]...))

    SnpData(people, snps, snparray, snp_info, person_info, des * ".bed", bimfile, famfile)
end


"""
    hcat(A...;des="tmp_hcat_" * string(hcat_counter))

Concatenate SnpData along dimension 2.
"""
hcat_counter = 1
function Base.hcat(A::SnpData...; des="tmp_hcat_" * string(hcat_counter))
    global hcat_counter
    if des == "tmp_hcat_" * string(hcat_counter)
        hcat_counter = hcat_counter + 1
    end

    # get person_info
    @assert allequal([x.person_info for x in A]) "person_info are different"
    person_info = A[1].person_info

    # vcat snp_info
    snp_info = vcat([x.snp_info for x in A]...)

    # hcat snparray
    snparray = hcat([x.snparray for x in A]...; des=des)
    
    people, snps = size(person_info, 1), size(snp_info, 1)

    bimfile = des * ".bim"
    famfile = des * ".fam"

    writedlm(bimfile, hcat([snp_info[!, k] for k in SNP_INFO_KEYS]...))
    writedlm(famfile, hcat([person_info[!, k] 
                                for k in PERSON_INFO_KEYS]...))

    SnpData(people, snps, snparray, snp_info, person_info, des * ".bed", bimfile, famfile)
end

"""
   hvcat(rows::Tuple{Vararg{Int}}, values...; des="tmp_hvcat_" * string(hvcat_counter))

Horizontal and vertical concatenation in one call. 
"""
hvcat_counter = 1
function Base.hvcat(rows::Tuple{Vararg{Int}}, A::SnpData...; des="tmp_hvcat" * string(hvcat_counter))
    global hvcat_counter
    if des == "tmp_hvcat_" * string(hvcat_counter)
        hvcat_counter = hvcat_counter + 1
    end

    num_block_rows = length(rows)
    

    # collect person_info
    a = 1

    person_info = A[1].person_info

    for i = 2:num_block_rows
        person_info = vcat(person_info, A[a].person_info)
        a += rows[i]
    end    
    # collect snp_info
    snp_info = vcat([x.snp_info for x in A[1:rows[1]]]...)

    # hvcat snparray
    snparray = hvcat(rows, [x.snparray for x in A]...; des=des)

    people, snps = size(person_info, 1), size(snp_info, 1)

    bimfile = des * ".bim"
    famfile = des * ".fam"

    writedlm(bimfile, hcat([snp_info[!, k] for k in SNP_INFO_KEYS]...))
    writedlm(famfile, hcat([person_info[!, k] 
                                for k in PERSON_INFO_KEYS]...))

    SnpData(people, snps, snparray, snp_info, person_info, des * ".bed", bimfile, famfile)
end

"""
    merge_plink(d)
    merge_plink(des, d)

Merge the SnpData of the splitted plink files. 
Returns: merged SnpData. If `des` is given, it is written on that destination.
"""
@inline function isless_chromosome(x::AbstractString, y::AbstractString)
    x_int = try parse(Int32, x) 
    catch 
        typemax(Int32)
    end
    y_int = try parse(Int32, y)
    catch
        typemax(Int32)
    end
    x_int == y_int ? x < y : x_int < y_int
end

function merge_plink(d::Dict{AbstractString, SnpData})
    ks = sort(collect(keys(d)), lt=isless_chromosome)
   
    # get person_info
    person_info = d[ks[1]].person_info
    
    # vcat snp_info
    @time snp_info = Base.vcat([d[k].snp_info for k in ks]...)
    # hcat snparray
    @time data = Base.hcat([d[k].snparray.data for k in ks]...)
    rowcounts = d[ks[1]].snparray.rowcounts
    @time columncounts = Base.hcat([d[k].snparray.columncounts for k in ks]...)
    snparray = SnpArray(data, rowcounts, columncounts, size(data)[1])

    people, snps = size(person_info,1), size(snp_info, 1)
    SnpData(people, snps, snparray, snp_info, person_info, "", "", "")
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
    @time for fn in matching_srcs
        chrsnpdata = SnpData(fn)
        chr = chrsnpdata.snp_info[!, :chromosome][1]
        @assert all(chr .== chrsnpdata.snp_info[!, :chromosome]) "Not all chrs are the same in $fn.bim."
        d[String(chr)] = chrsnpdata
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
    writedlm(bimfile, hcat([snp_info[!, k] for k in SNP_INFO_KEYS]...))
    
    # write fam file
    person_info = snpdata.person_info
    writedlm(famfile, hcat([person_info[!, k] 
                                for k in PERSON_INFO_KEYS]...))
    
    # write bed file
    open(bedfile, "w") do io
        write(io, 0x1b6c)
        write(io, 0x01)
        write(io, snpdata.snparray.data)
    end
    SnpData(snpdata.people, snpdata.snps, snpdata.snparray, snpdata.snp_info, snpdata.person_info, bimfile, bedfile, famfile)
end
