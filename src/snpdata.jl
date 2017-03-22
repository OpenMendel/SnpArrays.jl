#---------------------------------------------------------------------------#
# SnpData type
#---------------------------------------------------------------------------#

"""
Type for SNP and person information.
"""
type SnpData
  people::Int                       # number of rows (individuals)
  snps::Int                         # number of columns (snps)
  personid::Vector{AbstractString}  # names of individuals
  snpid::Vector{AbstractString}     # SNP ids
  chromosome::Vector{AbstractString}# SNP chromosome
  genetic_distance::Vector{Float64} # genetic distance
  basepairs::Vector{Int}            # SNP base pair position
  allele1::Vector{AbstractString}   # A1 code
  allele2::Vector{AbstractString}   # A2 code
  maf::Vector{Float64}              # minor allele frequencies
  minor_allele::BitVector           # bit vector designating the minor allele
  snpmatrix::SnpLike{2}             # matrix of genotypes or haplotypes
  missings_per_person::Vector{Int}  # number of missing genotypes per person
  missings_per_snp::Vector{Int}     # number of missing genotypes per snp
end # end SnpData

"""
Construct a SnpData type from a PLINK file.
"""
function SnpData(plink_file::AbstractString)

  # load snp info
  plink_bim_file = string(plink_file, ".bim")
  snp_info = readdlm(plink_bim_file, AbstractString)
  chromosome = snp_info[:, 1]
  snpid = snp_info[:, 2]
  genetic_distance = map(x -> parse(Float64, x), snp_info[:, 3])
  basepairs = map(x -> parse(Int, x), snp_info[:, 4])
  allele1 = snp_info[:, 5]
  allele2 = snp_info[:, 6]

  # load person info
  plink_fam_file = string(plink_file, ".fam")
  person_info = readdlm(plink_fam_file, AbstractString)
  personid = person_info[:, 2]

  # load snp array matrix
  snpmatrix = SnpArray(plink_file)
  maf, minor_allele, missings_per_snp, missings_per_person = summarize(snpmatrix)
  people, snps = size(snpmatrix)

  # construct SnpData unfiltered
  SnpData(people, snps, personid, snpid, chromosome, genetic_distance, basepairs,
    allele1, allele2, maf, minor_allele, snpmatrix, missings_per_person,
    missings_per_snp)
end

"""
Write snp data to Plink bed and bim files.
"""
function writeplink(filename::AbstractString, snpdata::SnpData)
  bimfile = filename * ".bim"
  bedfile = filename * ".bed"
  isfile(bimfile) && error("($bimfile) alread exists.")
  isfile(bedfile) && error("($bedfile) alread exists.")
  # write bim file
  writedlm(bimfile, zip(snpdata.chromosome, snpdata.snpid,
    snpdata.genetic_distance, snpdata.basepairs, snpdata.allele1, snpdata.allele2))
  # write bed file
  fid = open(bedfile, "w+")
  write(fid, UInt8[0x6c])
  write(fid, UInt8[0x1b])
  write(fid, UInt8[0x01])
  plinkbits = Mmap.mmap(fid, BitArray{3},
    (2, 4ceil(Int, 0.25snpdata.people), snpdata.snps))
  copy!(view(plinkbits, 1, 1:snpdata.people, :), snpdata.snpmatrix.A1)
  copy!(view(plinkbits, 2, 1:snpdata.people, :), snpdata.snpmatrix.A2)
  close(fid)
end

"""
    filter(snpdata, snp_idx, ppl_idx)

Filter a SnpData `snpdata` according SNP (column) index vector `snp_idx` and
person (column) index vector `ppl_idx`.

# Input
- `snpdata`: a SnpData.
- `snp_idx`: a Bitvector of SNP index.
- `ppl_idx`: a BitVector of person index.

# Output
- filtered snpdata.
"""
function filter(
  snpdata :: SnpData,
  snp_idx :: BitVector,
  ppl_idx :: BitVector
  )
  # subset vectors and matrices
  snpmatrix        = snpdata.snpmatrix[ppl_idx, snp_idx]
  personid         = snpdata.personid[ppl_idx]
  snpid            = snpdata.snpid[snp_idx]
  chromosome       = snpdata.chromosome[snp_idx]
  genetic_distance = snpdata.genetic_distance[snp_idx]
  basepairs        = snpdata.basepairs[snp_idx]
  allele1          = snpdata.allele1[snp_idx]
  allele2          = snpdata.allele2[snp_idx]

  # compute summary fields with summarize
  maf, minor_allele, missings_per_snp, missings_per_person = summarize(snpmatrix)

  # compute size fields with size
  people, snps = size(snpmatrix)

  SnpData(people, snps, personid, snpid, chromosome, genetic_distance,
    basepairs, allele1, allele2, maf, minor_allele, snpmatrix,
    missings_per_person, missings_per_snp)
end
