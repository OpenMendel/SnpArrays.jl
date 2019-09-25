"""
    reorder!(s::SnpData, order::Vector{<:Integer}; famnm=s.srcfam * ".reordered")
reorders the SnpData in-place, including the saved files if the sources exist. overwrites the fam file if `writefam` is true. 
SnpData should be opened with "r+" mode (e.g. `Snpdata("mouse", "r+")`).
WARNING: the data is directly modified on the disk.
"""
function reorder!(s::SnpData, order::Vector{<:Integer}; famnm=split(s.srcfam, ".fam")[1] * ".reordered.fam")
    s.snparray[:,:] .= s.snparray[order, :]
    s.person_info[:, :] .= s.person_info[order, :]
    writedlm(famnm, hcat([s.person_info[!, k] 
                            for k in PERSON_INFO_KEYS]...))
end