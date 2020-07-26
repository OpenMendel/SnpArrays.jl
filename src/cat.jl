@inline function allequal(x)
    length(x) < 2 && return true
    first_el = x[1]
    i = 2
    @inbounds for i=2:length(x)
        x[i] == first_el || return false
    end
    return true
end

"""
    vcat(A::SnpArray...;des="tmp_vcat_" * string(vcat_counter))

Concatenate SnpArray along dimension 1. (i.e. stack persons)
Only creates ".bed" file. Use vcat(A::SnpData...) for ".fam" and ".bim" files.
"""

arr_vcat_counter = 1
function Base.vcat(A::SnpArray...; des="tmp_vcat_arr_" * string(arr_vcat_counter))
    global arr_vcat_counter
    if des == "tmp_vcat_arr_" * string(arr_vcat_counter)
        arr_vcat_counter = arr_vcat_counter + 1
    end
    num_rows = [size(x)[1] for x in A]
    num_cols = [size(x)[2] for x in A]
    
    @assert allequal(num_cols) "number of columns are not the same"
    des_rows = sum(num_rows)
    des_cols = num_cols[1]
    cum_rows = Base.vcat([0], cumsum(num_rows))

    desbedfile = des * ".bed"

    makestream(desbedfile, "w+") do io
        write(io, 0x1b6c)
        write(io, 0x01)
        write(io, Matrix{UInt8}(undef, (des_rows + 3) >> 2, des_cols))
    end
    des_arr = SnpArray(desbedfile, des_rows, "r+")

    for i in 1:length(A)
        des_arr[(cum_rows[i] + 1): cum_rows[i + 1], :] .= @view A[i][:, :]
    end
    des_arr
end


"""
    hcat(A...;des="tmp_hcat_" * string(hcat_counter))

Concatenate SnpData along dimension 2. (i.e. stack locations)
Only creates ".bed" file. Use hcat(A::SnpData...) for ".fam" and ".bim" files.
"""
arr_hcat_counter = 1
function Base.hcat(A::SnpArray...; des="tmp_hcat_arr_" * string(arr_hcat_counter))
    global arr_hcat_counter
    if des == "tmp_hcat_arr_"*string(arr_vcat_counter)
        arr_hcat_counter = arr_hcat_counter + 1
    end
    num_rows = [size(x)[1] for x in A]
    num_cols = [size(x)[2] for x in A]

    @assert allequal(num_rows) "number of rows are not the same"
    des_rows = num_rows[1]
    des_cols = sum(num_cols)
    cum_cols = Base.vcat([0], cumsum(num_cols))

    desbedfile = des * ".bed"

    makestream(desbedfile, "w+") do io
        write(io, 0x1b6c)
        write(io, 0x01)
        write(io, Matrix{UInt8}(undef, (des_rows + 3) >> 2, des_cols))
    end
    des_arr = SnpArray(desbedfile, des_rows, "r+")

    des_arr.data .= hcat([s.data for s in A]...)
    des_arr
end

"""
   hvcat(rows::Tuple{Vararg{Int}}, values...; des="tmp_hvcat_arr_" * string(arr_hvcat_counter))

Horizontal and vertical concatenation in one call. 
Only creates ".bed" file. Use hvcat(rows::Tuple{Vararg{Int}}, values::SnpData...) for ".fam" and ".bim" files.
"""
arr_hvcat_counter = 1
function Base.hvcat(rows::Tuple{Vararg{Int}}, A::SnpArray...; des="tmp_hvcat_arr_" * string(arr_hvcat_counter))
    global arr_hvcat_counter
    if des == "tmp_hvcat_arr_"*string(arr_hvcat_counter)
        arr_hvcat_counter = arr_hvcat_counter + 1
    end

    num_block_rows = length(rows)

    num_cols = 0
    for i = 1:rows[1]
        num_cols += size(A[i], 2)
    end

    num_rows = 0
    a = 1
    for i = 1:num_block_rows
        num_rows += size(A[a], 1)
        a += rows[i]
    end

    desbedfile = des * ".bed"

    makestream(desbedfile, "w+") do io
        write(io, 0x1b6c)
        write(io, 0x01)
        write(io, Matrix{UInt8}(undef, (num_rows +3) >> 2, num_cols))
    end
    des_arr = SnpArray(desbedfile, num_rows, "r+")

    a = 1
    r = 1
    for i = 1:num_block_rows
        c = 1
        szi = size(A[a], 1)
        for j = 1:rows[i]
            Aj = A[a+j-1]
            szj = size(Aj, 2)
            if size(Aj, 1) != szi
                throw(ArgumentError("mismatched height in block row $(i) (expected $szi, got $(size(Aj, 1)))"))
            end
            if c-1+szj > num_cols
                throw(ArgumentError("block row $(i) has mismatched number of columns (expected $num_cols, got $(c-1+szj)"))
            end
            des_arr[r:r-1+szi, c:c-1+szj] .= @view Aj[:,:]
            c += szj
        end
        if c != num_cols+1
            throw(ArgumentError("block row $(i) has mismatched number of columns (expected $num_cols, got $(c-1))"))
        end
        r += szi
        a += rows[i]
    end
    des_arr 
end

