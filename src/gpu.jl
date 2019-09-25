using .OpenCL
export SnpVariables, SnpCLVariables

const Float = Union{Float64, Float32}
const gpucode64 = joinpath(@__DIR__, "kernels", "kernels64.cl")
const gpucode32 = joinpath(@__DIR__, "kernels", "kernels32.cl")
const cldict = Dict(Float64 => gpucode64, Float32 => gpucode32)


struct SnpVariables{T}
    model::Union{Val{1}, Val{2}, Val{3}}
    center::Bool
    scale::Bool
    μ::Vector{T}
    σinv::Vector{T}
    storagev1::Vector{T}
    storagev2::Vector{T}
end  

function SnpVariables{T}(s::AbstractSnpArray; model = ADDITIVE_MODEL, 
        center::Bool=false, scale::Bool = false) where T <: AbstractFloat
    if model == ADDITIVE_MODEL
        if center || scale
            μ = Vector{T}(undef, size(s, 2))
            μ[:] = mean(s, dims=1, model=ADDITIVE_MODEL)
        else
            μ = T[]
        end
        if scale
            σinv = Vector{T}(undef, size(s, 2))
            @inbounds @simd for j in 1:size(s, 2)
                σinv[j] = sqrt(μ[j] * (1 - μ[j] / 2))
                σinv[j] = σinv[j] > 0 ? inv(σinv[j]) : one(T)
            end
        else
            σinv = T[]
        end
    elseif model == DOMINANT_MODEL
        if center || scale
            μ = Vector{T}(undef, size(s, 2))
            μ[:] = mean(s, dims=1, model=DOMINANT_MODEL)
        else
            μ = T[]
        end
        if scale
            σinv = Vector{T}(undef, size(s, 2))
            @inbounds @simd for j in 1:size(s, 2)
                σinv[j] = sqrt(μ[j] * (1 - μ[j]))
                σinv[j] = σinv[j] > 0 ? inv(σinv[j]) : one(T)
            end
        else
            σinv = T[]
        end
    elseif model == RECESSIVE_MODEL
        if center || scale
            μ = Vector{T}(undef, size(s, 2))
            μ[:] = mean(s, dims=1, model=RECESSIVE_MODEL)
        else
            μ = T[]
        end
        if scale
            σinv = Vector{T}(undef, size(s, 2))
            @inbounds @simd for j in 1:size(s, 2)
                σinv[j] = sqrt(μ[j] * (1 - μ[j]))
                σinv[j] = σinv[j] > 0 ? inv(σinv[j]) : one(T)
            end
        else
            σinv = T[]
        end
    else
        throw(ArgumentError("unrecognized model $model"))
    end
    storagev1 = Vector{T}(undef, size(s, 1))
    storagev2 = Vector{T}(undef, size(s, 2))
    SnpVariables{T}(model, center, scale, μ, σinv, storagev1, storagev2)
end

struct SnpCLVariables{T <: Float, V <: cl.Buffer}
    df_buff     :: V
    x_buff      :: cl.Buffer{UInt8}
    y_buff      :: V
    mask_buff   :: cl.Buffer{Int}
    device      :: cl.Device
    ctx         :: cl.Context
    queue       :: cl.CmdQueue
    μ_buff      :: V
    σinv_buff      :: V
    red_buff    :: V
    xtyk        :: cl.Kernel
    rxtyk       :: cl.Kernel
    reset_x     :: cl.Kernel
    wg_size     :: Int
    y_chunks    :: Int
    r_chunks    :: Int
    m           :: Int
    n           :: Int
    #p2          :: Int 
    m32         :: Int32
    n32         :: Int32
    y_chunks32  :: Int32
    blocksize32 :: Int32
    wg_size32   :: Int32
    y_blocks32  :: Int32
    r_length32  :: Int32
    genofloat   :: cl.LocalMem{T}
end

function SnpCLVariables(
    df_buff     :: cl.Buffer{T},
    x_buff      :: cl.Buffer{UInt8},
    y_buff      :: cl.Buffer{T},
    mask_buff   :: cl.Buffer{Int},
    device      :: cl.Device,
    ctx         :: cl.Context,
    queue       :: cl.CmdQueue,
    μ_buff      :: cl.Buffer{T},
    σinv_buff      :: cl.Buffer{T},
    red_buff    :: cl.Buffer{T},
    xtyk        :: cl.Kernel,
    rxtyk       :: cl.Kernel,
    reset_x     :: cl.Kernel,
    wg_size     :: Int,
    y_chunks    :: Int,
    r_chunks    :: Int,
    m           :: Int,
    n           :: Int,
    #p2          :: Int,
    m32         :: Int32,
    n32         :: Int32,
    y_chunks32  :: Int32,
    blocksize32 :: Int32,
    wg_size32   :: Int32,
    y_blocks32  :: Int32,
    r_length32  :: Int32,
    genofloat   :: cl.LocalMem{T}
) where {T <: Float}
    
    #SnpCLVariables{T, eltype(y_buff)}(df_buff, x_buff, y_buff, mask_buff, device, ctx, queue, m_buff, p_buff, red_buff, xtyk, rxtyk, reset_x, wg_size, y_chunks, r_chunks, n, p, p2, n32, p32, y_chunks32, blocksize32, wg_size32, y_blocks32, r_length32, genofloat)
    SnpCLVariables{T, eltype(y_buff)}(df_buff, x_buff, y_buff, mask_buff, device, ctx, queue, μ_buff, σinv_buff, red_buff, xtyk, rxtyk, reset_x, wg_size, y_chunks, r_chunks, m, n, m32, n32, y_chunks32, blocksize32, wg_size32, y_blocks32, r_length32, genofloat)
end


function SnpCLVariables(
    z      :: DenseVector{T},
    x      :: AbstractSnpArray,
    y      :: DenseVector{T},
    #μ      :: DenseVector{T},
    #σinv   :: DenseVector{T},
    #x_bm :: SnpBitMatrix{T} = SnpBitMatrix{T}(x; model=ADDITIVE_MODEL, center=true, scale=true),
    kern   :: String = cldict[T],
    rmask :: Vector{Int} = ones(Int, size(y));
    model = ADDITIVE_MODEL,
    center = true,
    scale = true,
    xvars::SnpVariables{T} = SnpVariables{T}(x; model = model, center = center, scale = scale)
) where {T <: Float}
    #n,p         = size(x.geno)
    
    !(center && scale && model==ADDITIVE_MODEL) && throw(ArgumentError("Unsupported option: only centered, scaled, additive model is supported for now."))
    m, n          = size(x)
    
    wg_size     = 256
    y_chunks    = div(m, wg_size) + (m % wg_size != 0 ? 1 : 0)
    y_blocks    = div(y_chunks, wg_size) + (y_chunks % wg_size != 0 ? 1 : 0)
    r_chunks    = div(n*y_chunks, wg_size) + ((n*y_chunks) % wg_size != 0 ? 1 : 0)
    wg_size32   = convert(Int32, wg_size)
    m32         = convert(Int32, m)
    n32         = convert(Int32, n)
    y_chunks32  = convert(Int32, y_chunks)
    y_blocks32  = convert(Int32, y_blocks)
    # blocksize32 = convert(Int32, x.geno.blocksize)
    blocksize32 = convert(Int32, (x.m + 3) >> 2)
    r_length32  = convert(Int32, n*y_chunks)
    device      = last(cl.devices(:gpu))
    ctx         = cl.Context(device) :: cl.Context
    queue       = cl.CmdQueue(ctx) :: cl.CmdQueue
    program     = cl.Program(ctx, source=read(kern)) :: cl.Program
    cl.build!(program)
    xtyk        = cl.Kernel(program, "compute_xt_times_vector")
    rxtyk       = cl.Kernel(program, "reduce_xt_vec_chunks")
    reset_x     = cl.Kernel(program, "reset_x")
    #x_buff      = cl.Buffer(Int8, ctx, (:r,  :copy), hostbuf = sdata(x.geno.x)) :: cl.Buffer{Int8}
    x_buff      = cl.Buffer(UInt8, ctx, (:r, :copy), hostbuf = x.data) :: cl.Buffer{UInt8}
    y_buff      = cl.Buffer(T,    ctx, (:r,  :copy), hostbuf = y) :: cl.Buffer{T}
    #m_buff      = cl.Buffer(T,    ctx, (:r,  :copy), hostbuf = sdata(x.means)) :: cl.Buffer{T}
    μ_buff      = cl.Buffer(T,    ctx, (:r,  :copy), hostbuf = xvars.μ) :: cl.Buffer{T}
    #p_buff      = cl.Buffer(T,    ctx, (:r,  :copy), hostbuf = sdata(x.precs)) :: cl.Buffer{T}
    σinv_buff      = cl.Buffer(T,    ctx, (:r,  :copy), hostbuf = xvars.σinv) :: cl.Buffer{T}
    df_buff     = cl.Buffer(T,    ctx, (:rw, :copy), hostbuf = z) :: cl.Buffer{T}
    red_buff    = cl.Buffer(T,    ctx, (:rw), n*y_chunks) :: cl.Buffer{T}
    mask_buff   = cl.Buffer(Int,  ctx, (:r,  :copy), hostbuf = rmask) :: cl.Buffer{Int}
    genofloat   = cl.LocalMem(T, wg_size)

    #SnpCLVariables{T, cl.Buffer{T}}(df_buff, x_buff, y_buff, mask_buff, device, ctx, queue, m_buff, p_buff, red_buff, xtyk, rxtyk, reset_x, wg_size, y_chunks, r_chunks, n, p, 0, n32, p32, y_chunks32, blocksize32, wg_size32, y_blocks32, r_length32, genofloat)
    SnpCLVariables{T, cl.Buffer{T}}(df_buff, x_buff, y_buff, mask_buff, device, ctx, queue, μ_buff, σinv_buff, red_buff, xtyk, rxtyk, reset_x, wg_size, y_chunks, r_chunks, m, n, m32, n32, y_chunks32, blocksize32, wg_size32, y_blocks32, r_length32, genofloat)

end

function copy_y!(v::SnpCLVariables{T}, y::DenseVector{T}) where T <: Float
    cl.copy!(v.queue, v.y_buff, y) :: cl.NannyEvent
end

function reset_x!(v::SnpCLVariables{T}) where T <: Float
    v.queue(v.reset_x, (v.wg_size*v.r_chunks,1,1), nothing, v.red_buff, v.r_length32, zero(T)) :: cl.Event
end

function xty!(v::SnpCLVariables{T}) where T <: Float
    #cl.call(v.queue, v.xtyk, (v.wg_size*v.y_chunks, v.p, 1), (v.wg_size, 1, 1), v.n32, v.p32, v.y_chunks32, v.blocksize32, v.wg_size32, v.x_buff, v.red_buff, v.y_buff, v.m_buff, v.p_buff, v.mask_buff, v.genofloat) :: cl.Event
    v.queue(v.xtyk, (v.wg_size*v.y_chunks, v.n, 1), (v.wg_size, 1, 1), v.m32, v.n32, v.y_chunks32, v.blocksize32, v.wg_size32, v.x_buff, v.red_buff, v.y_buff, v.μ_buff, v.σinv_buff, v.mask_buff, v.genofloat) :: cl.Event
end

function xty_reduce!(v::SnpCLVariables{T}) where T <: Float
    #cl.call(v.queue, v.rxtyk, (v.wg_size,v.p,1), (v.wg_size,1,1), v.n32, v.y_chunks32, v.y_blocks32, v.wg_size32, v.red_buff, v.df_buff, v.genofloat) :: cl.Event
    v.queue(v.rxtyk, (v.wg_size,v.n,1), (v.wg_size,1,1), v.m32, v.y_chunks32, v.y_blocks32, v.wg_size32, v.red_buff, v.df_buff, v.genofloat) :: cl.Event
end

function copy_xty!(xty::DenseVector{T}, v::SnpCLVariables{T}) where T <: Float
    cl.copy!(v.queue, xty, v.df_buff) :: cl.NannyEvent
    #cl.queue(cl.copy!, sdata(xty), v.df_buff) :: cl.NannyEvent
end


"""
    LinearAlgebra.mul!(Xty::AbstractVector{T}, Xt::Transpose{UInt8, SnpArray}, 
        y::AbstractVector{T}; v::SnpClVariables)

Multiplies Xt and y using an OpenCL device.
"""
function LinearAlgebra.mul!(Xty :: AbstractVector{T}, Xt :: Transpose{UInt8, SnpArray}, y :: AbstractVector{T}; 
        v :: SnpCLVariables = SnpCLVariables(Xty, transpose(Xt), y),
        rmask :: Vector{Int} = ones(Int, size(y))) where T <: Float
    copy_y!(v, y)
    reset_x!(v)
    xty!(v)
    xty_reduce!(v)
    copy_xty!(Xty, v)
    return Xty
end

