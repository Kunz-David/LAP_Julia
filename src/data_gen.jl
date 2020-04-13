module data_gen

using LAP_julia
"""
    gen_rand_flow(flow_size, max_magnitude, tile_size, filter_amp)

gen_rand_flow generates random flow using a uniform distribution

# Inputs:
    - `flow_size`       -> dimensions of the flow
    - `max_magnitude`   -> maximum amplitude of the displacement
    - `tile_size`       -> size random element of start matrix (the larger the slower the flow)
    - `filter_amp`      -> size of the gaussian filter which is used to smooth the random start matrix
# Output:
    - `rand_flow`         -> generated random flow
"""
function gen_rand_flow(flow_size, max_magnitude, tile_size=Nothing, filter_amp=Nothing)

    # set default values
    if tile_size == Nothing
        tile_size = ceil(Int64, flow_size[1]/6)
    end
    if filter_amp == Nothing
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
    rand_flow = LAP_julia.helpers.clean_using_gaussain(rand_flow, [filter_amp, filter_amp])

    return rand_flow

end


end # module
