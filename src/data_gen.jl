module data_gen

export gen_rand_flow, gen_chess, gen_init

using TestImages
using LAP_julia


"""
    gen_rand_flow(flow_size::Tuple{T, T}=(200, 200), max_magnitude::Real=20, tile_size=nothing; filter_amp=nothing)::Flow where {T <: Integer}

Generate a smoothly varying random flow. The flow parameters are set by the function arguments.

It works by generating a tiled flow where each tile has a random uniform flow. Then it smooths these tiles by filtering with a gaussian.

# Arguments
- `flow_size::Tuple{T, T}=(200, 200)`: dimensions of the flow.
- `max_magnitude::Real=20`: maximum allowed amplitude of the displacement.
- `tile_size=nothing`: size of random uniform flow tiles that make up the start matrix. (The larger the slower the flow.) Note: If set to `flow_size` or larger it will generate a uniform pixel shift in a random direction.
- `filter_amp=nothing`: size of the gaussian filter which is used to smooth the random start matrix.

See also: [`showflow`](@ref), [`Flow`](@ref), [`gen_chess`](@ref)
"""
function gen_rand_flow(flow_size::Tuple{T, T}=(200, 200), max_magnitude::Real=20, tile_size=nothing, filter_amp=nothing)::Flow where {T <: Integer}

    # set default values
    if tile_size == nothing
        tile_size = ceil(Int64, flow_size[1]/6)
    end
    if filter_amp == nothing
        filter_amp = ceil(Int64, tile_size/2)
        if isodd(filter_amp)
            filter_amp += 1
        end
    end

    tile_count = ceil.(Int64, flow_size ./ tile_size)

    uv_base = Array{Float64}(undef, tile_count..., 2)
    uv_base[:,:,1] .= -max_magnitude .+ (2 .* max_magnitude .* rand(tile_count...))
    uv_base[:,:,2] .= -max_magnitude .+ (2 .* max_magnitude .* rand(tile_count...))

    uv_flow = Array{Float64}(undef, tile_size*tile_count[1], tile_size*tile_count[2], 2)
    uv_flow[:,:,1] .= repeat(uv_base[:,:,1], inner=(tile_size, tile_size))
    uv_flow[:,:,2] .= repeat(uv_base[:,:,2], inner=(tile_size, tile_size))

    # cut off what we dont want
    #rand_flow = zeros(Complex{Float64}, flow_size)
    rand_flow = uv_flow[1:flow_size[1], 1:flow_size[2], 1] .+ (im * uv_flow[1:flow_size[1], 1:flow_size[2], 2]);

    # blur to make it continuous
    rand_flow = LAP_julia.smooth_with_gaussian(rand_flow, [filter_amp, filter_amp])

    return rand_flow
end

"""
    gen_chess(tile_size::Integer=50, board_size::Integer=4)

Create a chessboard image with `board_size` tiles in each dimension, where each tile is `tile_size` pixels in each dimension.
Note that `board_size` has to be even.

See also: [`imgshow`](@ref), [`gen_rand_flow`](@ref)
"""
function gen_chess(tile_size::Integer=50, board_size::Integer=4)
    mini_board = [zeros(tile_size, tile_size) ones(tile_size, tile_size);
                  ones(tile_size, tile_size) zeros(tile_size, tile_size)]

    chessboard = repeat(mini_board, outer=(convert(Integer, (board_size/2)), convert(Integer, (board_size/2))))
    return chessboard
end

"""
    gen_init(type::Symbol=:lena; flow_args=[])

Create the usual testing data; img, imgw, flow
"""
function gen_init(type::Symbol=:lena; flow_args=[])
    if type == :lena
        img = testimage("lena_gray")
        img = Float64.(img)
    elseif type == :chess
        img = gen_chess(50,4)
    end

    flow = gen_rand_flow(size(img), flow_args...)
    imgw = warp_img(img, -real(flow), -imag(flow))

    return img, imgw, flow
end

end # module
