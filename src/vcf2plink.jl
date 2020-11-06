chrom_map = Dict{String, String}()
for i in 1:23
    chrom_map[string(i)] = string(i)
    chrom_map["chr" * string(i)] = string(i)
end
chrom_map["chrX"] = "X"
chrom_map["chrY"] = "Y"
chrom_map["chrXY"] = "XY"
chrom_map["chrMT"] = "MT"
chrom_map["chrM"] = "M"
chrom_map["X"] = "X"
chrom_map["Y"] = "Y"
chrom_map["XY"] = "XY"
chrom_map["MT"] = "MT"
chrom_map["M"] = "M"

@inline function geno_ismissing(record::VCF.Record, range::UnitRange{Int})
    return record.data[first(range)] == UInt8('.') || record.data[last(range)] == UInt8('.')
end

"""
    vcf2plink(vcffile, plinkprefix)

Convert VCF file to PLINK format. This drops multi-allelic variants and variants with missing ID. Reference alleles are set to A1.
"""
function vcf2plink(vcffile, plinkprefix)
    @assert endswith(vcffile, "vcf.gz") || endswith(vcffile, ".vcf")
    nsamples = begin
        vcfio = makestream(vcffile, "r")
        reader = VCF.Reader(vcfio)
        samples = length(VCF.header(reader).sampleID)
        makestream(plinkprefix * ".fam", "w") do io
            for (i, id) in enumerate(VCF.header(reader).sampleID)
                println(io, "$i\t$id\t0\t0\t0\t0")
            end
        end
        close(reader)
        samples
    end

    nsnps = begin
        vcfio = VCF.Reader(SnpArrays.makestream(vcffile, "r"))
        records = 0 
        for record in vcfio
            length(VCF.alt(record)) > 1 && continue
            try
                ismissing(VCF.id(record))
            catch e
                continue
            end
            records += 1
        end
        close(vcfio)
        records
    end

    makestream(plinkprefix * ".bed", "w+") do io
        write(io, 0x1b6c)
        write(io, 0x01)
        write(io, Matrix{UInt8}(undef, (nsamples + 3) >> 2, nsnps))
    end
    snparray = SnpArray(plinkprefix * ".bed", nsamples, "r+")
    bimio = makestream(plinkprefix * ".bim", "w")
    vcfio = makestream(vcffile, "r")
    reader = VCF.Reader(vcfio)

    j = 0
    for record in reader
        length(VCF.alt(record)) > 1 && continue
        try
            ismissing(VCF.id(record))
        catch e
            continue
        end
        j += 1
        chrom = VCF.chrom(record)
        id = VCF.id(record)[1]
        dist = 0
        pos = VCF.pos(record)
        ref = VCF.ref(record)
        alt = VCF.alt(record)[1]
        println(bimio, "$chrom\t$id\t$dist\t$pos\t$ref\t$alt")
        
        gtkey = VCF.findgenokey(record, "GT")
        for (i, _) in enumerate(VCF.header(reader).sampleID)
            geno = record.genotype[i]
            # dropped field or "." => 0x2e
            if gtkey > lastindex(geno) || geno_ismissing(record, geno[gtkey])
                snparray[i,j] = 0x01
            else
                # "0" => 0x30, "1" => 0x31
                if record.data[geno[gtkey][1]] == 0x30
                    if record.data[geno[gtkey][3]] == 0x30
                        snparray[i, j] = 0x00
                    else
                        snparray[i, j] = 0x02
                    end
                elseif record.data[geno[gtkey][1]] == 0x31
                    if record.data[geno[gtkey][3]] == 0x31
                        snparray[i, j] = 0x03
                    else
                        snparray[i, j] = 0x02
                    end
                end
            end
        end
    end
    close(bimio)
    close(vcfio)
    close(reader)
    snparray
end
