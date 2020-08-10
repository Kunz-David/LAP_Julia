
using ImageFiltering: centered, KernelFactors.gaussian, kernelfactors, imfilter
using ImageDistances
using TimerOutputs


"""
    move_landmarks(locations, flow)

Move each `location::Tuple(T,T) where T <: Real` from the `locations` vector using the `flow`.
"""
function move_landmarks(locations, flow)
    itp = Interpolations.interpolate(flow, BSpline(Linear()))
    shifted_locations = Array{Float64}(undef, size(locations))
    for k in 1:size(locations, 1)
        shifted_locations[k,:] .= locations[k,:] .+ [imag(itp(locations[k,:]...)), real(itp(locations[k,:]...))]
    end
    return shifted_locations
end

"""
    pad_images(image1::Image, image2::Image)

Adds zeros to the right and bottom of `image1` and `image2` to make them the same size.
"""
function pad_images(image1::Image, image2::Image)

    (a, b) = size(image1)
    (c, d) = size(image2)

    # if there is nothing to do, end
    if (a, b) == (c, d)
        return image1, image2
    end

    if (a < c)
        image1 = [image1; zeros(eltype(image1), c - a, b)];
        a = c
    elseif (a > c)
        image2 = [image2; zeros(eltype(image2), a - c, d)];
        c = a
    end

    if (b < d)
        image1 = [image1 zeros(eltype(image1), a, d - b)];
    elseif (b > d)
        image2 = [image2 zeros(eltype(image2), c, b - d)];
    end
    return image1, image2
end

"""
    pad_images(image1::Image, image2::Image)

Adds zeros to the right and bottom of `image1` and `image2` to make them the same size. Also create a mask,
which has zeros at the edges of added zeros.
"""
function pad_images_mask(target::Image, source::Image)

    (a, b) = size(target)
    (c, d) = size(source)

    mask = trues(a > c ? a : c, b > d ? b : d)

    # if there is nothing to do, end
    if (a, b) == (c, d)
        return target, source, mask
    end

    if (a < c)
        target = [target; zeros(eltype(target), c - a, b)];
        # add horizontal mask
        mask[a, :] .= false
        mask[a+1, :] .= false
        a = c
    elseif (a > c)
        source = [source; zeros(eltype(source), a - c, d)];
        c = a
    end

    if (b < d)
        target = [target zeros(eltype(target), a, d - b)];
        # add vertical mask
        mask[:, b] .= false
        mask[:, b+1] .= false
    elseif (b > d)
        source = [source zeros(eltype(source), c, b - d)];
    end
    return target, source, mask
end

"""
    normalize_to_zero_one(arr)

Normalize `arr` so that it has a maximum of 1 and minimum of 0.
"""
function normalize_to_zero_one(arr)
    normalized = arr .- minimum(arr)
    normalized = normalized ./ maximum(normalized)
    return normalized
end


"""
    rescale!(image1::Image, image2::Image)

Rescale `image1` and `image2` intensities to span the whole `[0, 1]`.
"""
function rescale!(image1::Image, image2::Image)

    max_intensity = maximum([image1 image2])
    min_intensity = minimum([image1 image2])

    image1 = (image1 .- min_intensity)./(max_intensity - min_intensity);
    image2 = (image2 .- min_intensity)./(max_intensity - min_intensity);

    return image1, image2
end

"""
    smooth_with_gaussian!(A::Matrix{<:Number}, window_half_size_one_dim::Integer)

Clean the Matrix `A` by smoothing using a square 2D Gaussian filter of size `2 * window_half_size_one_dim + 1` in each dimension.
"""
@inline function smooth_with_gaussian!(A::Matrix{<:Number}, window_half_size_one_dim::Integer)# where T<:Integer
    window_half_size = [window_half_size_one_dim, window_half_size_one_dim]
    return smooth_with_gaussian!(A, window_half_size)
end

"""
    smooth_with_gaussian!(A::Matrix{<:Number}, window_half_size)

Smooth the Matrix `A` by smoothing using a 2D Gaussian filter of size `2 * window_half_size + 1`.
"""
function smooth_with_gaussian!(A::Matrix{<:Number}, window_half_size; timer=TimerOutput("smoothing"))

    σ = 2 * window_half_size[1]

    @timeit_debug timer "create 1d gauss" begin
        gaus = gaussian(σ, 2 * σ + 1)
    end
    @timeit_debug timer "imfilter" begin
        imfilter!(A, A, kernelfactors((gaus, gaus)), "symmetric")
    end
    return A
end


function create_sparse_flow_from_sparse(flow_vectors, inds, flow_size)
    flow = zeros(flow_size) .+ zeros(flow_size) .* im
    for (vec, ind) in zip(flow_vectors, inds)
        flow[ind] .= vec
    end
    return flow
end

function create_sparse_flow_from_full(full_flow, inds)
    sparse_flow = zeros(size(full_flow)) .+ zeros(size(full_flow)) .* im

    for ind in inds
        sparse_flow[ind] .= full_flow[ind]
    end
    return sparse_flow
end


"""
    function angle_rmse(x, y)

Calculate the root mean square error in angle between `x` and `y`. Output in degrees.
"""
function angle_rmse(x, y)
    # @assert eltype(x) == eltype(y)
    @assert eltype(x) <: Complex

    return sqrt(mse(rad2deg.(angle.(x)), rad2deg.(angle.(y))))
end

"""
    function angle_mae(x, y)

Calculate the mean absolute error in angle between `x` and `y`. Output in degrees.
"""
function angle_mae(x, y)
    return mean(abs.(rad2deg.(angle.(x)) - rad2deg.(angle.(y))))
end

"""
    function vec_len(x)

Calculate the lenght of vector `x`. `x` is a complex number.
"""
function vec_len(x)
    return sqrt(real(x)^2 + imag(x)^2)
end

"""
    function mean(x)

Calculate the mean of `x`.
"""
function mean(x)
    sum(x)/length(x)
end

# TODO: sppeeudp
"""
    inds_to_points(inds::Array{CartesianIndex{2}, 1})

Transform an array of `CartersianIndexes` to an array of where each column is a vector of indices of the input array.
"""
function inds_to_points(inds::Array{CartesianIndex{2}, 1})
    pos_x = [ind[1] for ind in inds]
    pos_y = [ind[2] for ind in inds]
    return collect(transpose([pos_x pos_y]))
end


"""
    points_to_inds(points::Array{Int64,2})

Reverse of the `inds_to_points` function.
"""
function points_to_inds(points)
    return [CartesianIndex(points[1, k], points[2, k]) for k in 1:size(points,2)]
end

"""
    max_displacement(flow)

Find the maximum displacement of `flow`, ignoring `NaN`s.
"""
function max_displacement(flow)
    magnitudes = map(x -> vec_len(x), flow)
    max_mag = maximum(filter(!isnan, magnitudes))
    return max_mag
end


"""
    lap(args; kwargs)

Return a 2D `Flow` matrix of displacements that transforms `source` closer to `target` and `source` image transformed with this displacement field, `source_reg`.

# Details
Perform the classic Local All-Pass algorithm with post-proccessing (inpainting and smoothing), then warps the `source` image with the resulting estimate flow.

# Arguments:
- `target::Image`: target/fixed grayscale image.
- `source::Image`: source/moving grayscale image.
- `fhs::Integer`: the half size of the base of the gaussian filters used.
- `window_size=(2fhs+1, 2fhs+1)`: the size of the local window (tuple of 2 ints), usually same as filter_size.

# Keyword Arguments:
- `timer::TimerOutput=TimerOutput("lap")`: provide a timer which times certain blocks in the function.
- `display::Bool=false`: verbose and debug prints.

# Outputs:
- `flow_estim`: full estimated flow.
- `source_reg`: `source` image warped by `flow_estim`.

See also: [`single_lap`](@ref), [`inpaint_nans!`](@ref), [`smooth_with_gaussian!`](@ref), [`warp_img`](@ref).
"""
function lap(target::Image,
             source::Image,
             fhs,
             window_size=(2fhs+1, 2fhs+1);
             timer::TimerOutput=TimerOutput("lap"),
             display::Bool=false)

    @timeit_debug timer "lap" begin
        flow_estim = single_lap(target, source, fhs, window_size, timer=timer, display=display) # , new_feature=new_feature
    end
    @timeit_debug timer "inpainting" begin
        @timeit_debug timer "replicating borders" begin
            window_half_size = Int64.((window_size.-(1,1))./2)
            middle_vals = flow_estim[window_half_size[1]+1:end-window_half_size[1],
                                     window_half_size[2]+1:end-window_half_size[2]]
            flow_estim = parent(padarray(middle_vals, Pad(:replicate, window_half_size...)))
        end # "replication borders"
        @timeit_debug timer "inside" begin
            inpaint_nans!(flow_estim)
        end # "inside"
    end
    @timeit_debug timer "smoothing" begin
        smooth_with_gaussian!(flow_estim, window_half_size, timer=timer)
    end
    @timeit_debug timer "generate source_reg" begin
        source_reg = warp_img(source, -real(flow_estim), -imag(flow_estim))
        # source_reg = warp_img(source, -real(flow_estim), -imag(flow_estim), target)
    end
    # if display
    #     print_timer(timer)
    # end
    return flow_estim, source_reg
end

"""
    sparse_lap(args; kwargs)

Return a 2D `Flow` matrix of displacements that transforms `source` closer to `target`.

# Details
Uses the Local All-Pass algorithm at indices, `inds`, of the target image, with a high gradient magnitude, to estimate displacement vectors at these locations.
These estimated vectors are then fit into a parametric model (`flow_interpolation_method`) and a full displacemnt field/flow, or `full_flow_estim`, is made.
This flow is then used to warp the `source` image, creating the `source_reg` image.

# Arguments:
- `target::Image`: target/fixed grayscale image.
- `source::Image`: source/moving grayscale image.
- `fhs::Integer`: the half size of the base of the gaussian filters used.
- `window_size=(2fhs+1, 2fhs+1)`: the size of the local window (tuple of 2 ints), usually same as filter_size.

# Keyword Arguments:
- `spacing::Integer=25`: the smallest number of pixels that can separate two `inds`. See also: [`find_edge_points`](@ref)
- `point_count::Integer=100`: the number of `inds` we are looking for.
(The method will use the maximum amount up to this number if the this amount of `inds` isn't found.) See also: [`find_edge_points`](@ref)
- `timer::TimerOutput=TimerOutput("sparse_lap")`: provide a timer which times certain blocks in the function.
- `display::Bool=false`: verbose and debug prints.
- `flow_interpolation_method::Symbol=:quad`: choose which strategy to use for sparse flow interpolation.
(Choices: `:quad` -> fits to a global quadratic model, `:rbf` -> fits to a local rbf model.) See also: [`interpolate_flow`](@ref)
- `base_method_kwargs=Dict(:timer => timer, :display => display))`: keyword arguments passed to the base method.
(In this case: [`single_lap_at_points`](@ref)).

# Outputs:
- `full_flow_estim`: full estimated flow.
- `source_reg`: `source` image warped by `full_flow_estim`. See also: [`warp_img`](@ref).
- `flow_estim_at_inds`: estimated displacement vectors, at the indices `inds`. These are fitted into the `flow_interpolation_method` to
make the full estimated flow `full_flow_estim`.
- `inds`: indices chosen for the high gradient at which the `flow_estim_at_inds` is estimated.

See also: [`showflow`](@ref), [`imgshow`](@ref), [`find_edge_points`](@ref), [`interpolate_flow`](@ref)
"""
function sparse_lap(target,
                    source,
                    fhs,
                    window_size=(2fhs+1, 2fhs+1);
                    spacing::Integer=15,
                    point_count::Integer=200,
                    timer::TimerOutput=TimerOutput("sparse_lap"),
                    display::Bool=false,
                    flow_interpolation_method::Symbol=:quad,
                    base_method_kwargs=Dict(:timer => timer, :display => display))

    mask = parent(padarray(trues(size(target).-(2*fhs, 2*fhs)), Fill(false, (fhs, fhs), (fhs, fhs))))
    @timeit_debug timer "find edge points" begin
        inds = find_edge_points(target, spacing=spacing, point_count=point_count, mask=mask)
    end
    @timeit_debug timer "sparse lap" begin
        flow_estim_at_inds, inds = single_lap_at_points(target, source, fhs, window_size, inds; base_method_kwargs...)
    end
    if all(isnan, flow_estim_at_inds)
        @timeit_debug timer "interpolate flow" begin
            full_flow_estim = zeros(size(target)) .* im .+ zeros(size(target))
        end
    else
        @timeit_debug timer "interpolate flow" begin
            full_flow_estim = interpolate_flow(flow_estim_at_inds, inds, size(target), method=flow_interpolation_method)
        end
    end
    @timeit_debug timer "generate source_reg" begin
        # source_reg = warp_img(source, -real(full_flow_estim), -imag(full_flow_estim))
        source_reg = warp_img(source, -real(full_flow_estim), -imag(full_flow_estim), target)
    end
    return full_flow_estim, source_reg, flow_estim_at_inds, inds
end
