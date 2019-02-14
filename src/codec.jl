const ALLOWED_FORMAT = ["gz", "zlib", "zz", "xz", "zst", "bz2"]

"""
    makestream(filepath)

Open a file with an appropriate codec determined by file suffix.
"""
function makestream(filepath, args...; kwargs...)
    io = open(filepath, args...; kwargs...)
    if endswith(filepath, ".gz")
        codec = isreadonly(io) ? GzipDecompressor() : GzipCompressor()
    elseif endswith(filepath, ".zlib")
        codec = isreadonly(io) ? ZlibDecompressor() : ZlibCompressor()
    elseif endswith(filepath, ".zz")
        codec = isreadonly(io) ? DeflateDecompressor() : DeflateCompressor()
    elseif endswith(filepath, ".xz")
        codec = isreadonly(io) ? XzDecompressor() : XzCompressor()
    elseif endswith(filepath, ".zst")
        codec = isreadonly(io) ? ZstdDecompressor() : ZstdCompressor()
    elseif endswith(filepath, ".bz2")
        codec = isreadonly(io) ? Bzip2Decompressor() : Bzip2Compressor()
    else
        return io
    end
    TranscodingStream(codec, io)
end

function makestream(f::Function, args...)
    io = makestream(args...)
    try
        f(io)
    finally
        close(io)
    end
end

"""
    checkplinkfilename(filename, suffix)

Check a filename is valid plink file name. Suffix should be "bed", "fam" or "bim".    
"""
function checkplinkfilename(filename::AbstractString, suffix::AbstractString)
    suffix ∈ ["bed", "fam", "bim"] || 
    throw(ArgumentError("suffix should be bed, fam or bim"))
    endswith(filename, "." * suffix) || 
    any(endswith.(filename, ".$suffix." .* ALLOWED_FORMAT)) || 
    throw(ArgumentError("compressed format should be one of $ALLOWED_FORMAT"))
end

"""
    compresse_plink(plink_file, format[="gz"], outfile[=plink_file])

Compress a set of Plink files to `gz` (default), `zlib` or `zz` format.
"""
function compress_plink(
    plinkfile::AbstractString, 
    format::AbstractString="gz", 
    outfile::AbstractString=plinkfile
    )
    format ∈ ALLOWED_FORMAT || throw(ArgumentError("compressed format should be one of $ALLOWED_FORMAT"))
    for suffix in [".bed", ".fam", ".bim"]
        fname = plinkfile * suffix
        if isfile(fname)
            open(fname) do input
                makestream(fname * "." * format, "w") do output
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
    format ∈ ALLOWED_FORMAT || throw(ArgumentError("compressed format should be one of $ALLOWED_FORMAT"))
    for suffix in [".bed", ".fam", ".bim"]
        iname = plink_compressed * suffix * "." * format
        oname = outfile * suffix
        if isfile(iname)
            open(oname, "w") do output
                makestream(iname) do input
                    write(output, input)
                end
            end
        end
    end
end
