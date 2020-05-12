module lap

export polyfilter_lap, single_lap

using LAP_julia
using ImageFiltering: Fill, KernelFactors.gaussian, centered, kernelfactors, imfilter!, padarray, Pad
using LinearAlgebra: qr
using PyPlot: Figure

function test()
    print("testing revise")
end

# Function implements the standard Local All-Pass Filter (LAP) algorithm
# (i.e. no multi-scale). Given input images I1 and I2, and basis (the set
# of filters), the function calculates the optical flow between the two
# images. Currently designed to use grayscale images.
#
# Inputs:
#           I1      -> Input image 1, gray scale, s1 by s2 array
#           I2      -> Input image 2, gray scale, s1 by s2 array
#           basis   -> Elements of the filter basis. M by N array, where N
#                   is the number of basis elements and M = R^2 where is
#                   the filter is R by R in size.
#           K1     -> Size of filter basis
#           W      -> Size of the local window used in the algorithm.
#
# Outputs:
#           uest    -> Estimate of the optical flow as a complex number
#           coeffs  -> Estimated coefficients for each filter basis
#           err     -> Error for each pixel: ||p(-x)*I2(x) - p(x)*I1(x)||^2
#

"""
    single_lap(image_1::Image, image_2::Image, filter_half_size::Integer, window_size, filter_num::Integer=3)

Return an estimate of a smoothly varying flow of size of `image_1` which is the displacement that transforms `image_2` closer to `image_1`.

# Arguments:
- `image_1::Image`: grayscale image 1.
- `image_2::Image`: grayscale image 2.
- `filter_num::Integer=3`: the number of basis filters used (so far only =3 implemented).
- `filter_half_size`: the half size of the base of the gaussian filters used.
- `window_size`: the size of the local window (list of 2 ints) usually same as filter_size.

See also: [`polyfilter_lap`](@ref), [`showflow`](@ref), [`warp_imgshowflow`](@ref), [`imgshowflow`](@ref)
"""
function single_lap(image_1::Image, image_2::Image, filter_half_size::Integer, window_size, filter_num::Integer=3)

    pixel_count = length(image_1)
    image_size = size(image_1)
    filter_size = 2 * filter_half_size + 1

    # Prepare filters:
    # Calculate separable filters from basis:
    sigma = (filter_half_size + 2) / 4
    centered_inds = centered(-filter_half_size:filter_half_size)
    gaus = gaussian(sigma, filter_size)
    gaus_der = gaus .* centered_inds .* (-1)/sigma^2
    gaus_der_flip = reverse(gaus_der, dims=1)

    # prepare 2D filter basis for later use:
    basis = zeros(filter_size, filter_size, filter_num)

    # dif = (image forward - image backward) of each filter
    dif = zeros(image_size[1], image_size[2], filter_num)

    if filter_num == 3
        # temporary place to store filtered images
        tmp_filtered_1 = similar(image_1)
        tmp_filtered_2 = similar(image_2)

        # basis 1 - (gaus, gaus) both forward and backward
        kernf_1 = kernelfactors((gaus, gaus))
        basis[:, :, 1] = broadcast(*, kernf_1...)
        imfilter!(tmp_filtered_1, image_1, kernf_1, "symmetric")
        imfilter!(tmp_filtered_2, image_2, kernf_1, "symmetric")
        dif[:, :, 1] = tmp_filtered_2 - tmp_filtered_1

        # basis 2 - (gaus, gaus_der_flip) as forward, (gaus, gaus_der) as backward
        kernf_2f = kernelfactors((gaus, gaus_der_flip))
        basis[:, :, 2] = broadcast(*, kernf_2f...)
        imfilter!(tmp_filtered_1, image_1, kernf_2f, "symmetric")
        kernf_2b = kernelfactors((gaus, gaus_der))
        imfilter!(tmp_filtered_2, image_2, kernf_2b, "symmetric")
        dif[:, :, 2] = tmp_filtered_2 - tmp_filtered_1

        # basis 3 - (gaus_der_flip, gaus) as forward, (gaus_der, gaus) as backward
        kernf_3f = kernelfactors((gaus_der_flip, gaus))
        basis[:, :, 3] = broadcast(*, kernf_3f...)
        imfilter!(tmp_filtered_1, image_1, kernf_3f, "symmetric")
        kernf_3b = kernelfactors((gaus_der, gaus))
        imfilter!(tmp_filtered_2, image_2, kernf_3b, "symmetric")
        dif[:, :, 3] = tmp_filtered_2 - tmp_filtered_1
    end

    dif = reshape(dif, (:, filter_num))

    # prepare matrices for linear system of equations
    A = zeros(filter_num-1, filter_num-1, pixel_count)
    b = zeros(filter_num-1, pixel_count)

    for k in 1:filter_num-1
        for l in k:filter_num-1
            @views window_sum!(A[k, l, :], dif[:, k+1] .* dif[:, l+1], image_size, window_size)
            A[l, k, :] = A[k, l, :]
        end
        @views window_sum!(b[k, :], dif[:, k+1] .* dif[:, 1] .* (-1), image_size, window_size)
    end

    # Perform Gauss elimination on all pixels in parallel:
    # coeffs will be of shape: pixel_count, filter_num-1
    @views coeffs = multi_mat_div(A, b)
    # adding ones so that all base filters have their coefficients even the first one
    all_coeffs = [ones(pixel_count) coeffs]

    # make a border mask. Is 0 if its in a border of filter_half_size size.
    window_half_size = Int64.((window_size .- 1) ./ 2)
    border_mask = parent(padarray(ones((image_size .- (2 .* window_half_size))...),
                    Fill(NaN, window_half_size, window_half_size)))
    @views all_coeffs = all_coeffs .* reshape(border_mask, (:, 1))

    k = (-filter_half_size:filter_half_size)

    # Get the displacement vector field from the filters
    u1_top = zeros(image_size);
    u1_bot = zeros(image_size);
    u2_top = zeros(image_size);
    u2_bot = zeros(image_size);

    for n in 1:filter_num
        @views u1_top[:] = u1_top[:] .- sum(transpose(basis[:, :, n]) * k) .* all_coeffs[:, n];
        @views u1_bot[:] = u1_bot[:] .+ sum(basis[:, :, n]) .* all_coeffs[:, n];

        @views u2_top[:] = u2_top[:] .- sum(basis[:, :, n] * k) .* all_coeffs[:, n];
        @views u2_bot[:] = u2_bot[:] .+ sum(basis[:, :, n]) .* all_coeffs[:, n];
    end

    u_est = 2 .* ((im .* u1_top ./ u1_bot) .+ (u2_top ./ u2_bot));

    # dont use estimations whose displacement is larger than the filter_half_size
    displacement = real(u_est).^2 + imag(u_est).^2
    displacement_mask = displacement .> filter_half_size^2
    u_est[displacement_mask] .= NaN .+ NaN .* 1im;

    return u_est
end

using ComputationalResources, Images

# function filt_onebyone!(imgfilt, img, kernel, filter_half_size, points)
#     padded = padarray(img, Pad(:symmetric, filter_half_size, filter_half_size))
#
#     @simd for k in 1:size(points, 2)
#         ind_filt = (points[1, k]:points[1, k], points[2, k]:points[2, k])
#         imfilter!(CPU1(Algorithm.FIRTiled()), imgfilt, padded, kernel, NoPad(), ind_filt)
#     end
# end

function filt_onebyone!(imgfilt, img, kernel, fhs, points)
    padded = padarray(img, Pad(:symmetric, fhs, fhs))

    @simd for k in 1:size(points, 2)
        ind_filt = (points[1, k]-fhs:points[1, k]+fhs, points[2, k]-fhs:points[2, k]+fhs)
        imfilter!(CPU1(Algorithm.FIRTiled()), imgfilt, padded, kernel, NoPad(), ind_filt)
    end
end


function single_lap(image_1::Image, image_2::Image, filter_half_size::Integer, window_size, filter_num, points)

    image_size = size(image_1)
    pixel_count = length(image_1)
    point_count = size(points, 2)
    filter_size = 2 * filter_half_size + 1

    # Prepare filters:
    # Calculate separable filters from basis:
    sigma = (filter_half_size + 2) / 4
    centered_inds = centered(-filter_half_size:filter_half_size)
    gaus = gaussian(sigma, filter_size)
    gaus_der = gaus .* centered_inds .* (-1)/sigma^2
    gaus_der_flip = reverse(gaus_der, dims=1)

    # prepare 2D filter basis for later use
    basis = zeros(filter_size, filter_size, filter_num)

    # dif = (points of forward image - points of backward image) for each filter
    dif = zeros(image_size..., filter_num)

    if filter_num == 3
        # temporary place to store filtered images
        tmp_filtered_1 = zeros(image_size)
        tmp_filtered_2 = zeros(image_size)

        # basis 1 - (gaus, gaus) both forward and backward
        kernf_1 = kernelfactors((gaus, gaus))
        basis[:, :, 1] = broadcast(*, kernf_1...)
        filt_onebyone!(tmp_filtered_1, image_1, kernf_1, filter_half_size, points)
        filt_onebyone!(tmp_filtered_2, image_2, kernf_1, filter_half_size, points)
        dif[:, :, 1] = tmp_filtered_2 - tmp_filtered_1

        # basis 2 - (gaus, gaus_der_flip) as forward, (gaus, gaus_der) as backward
        kernf_2f = kernelfactors((gaus, gaus_der_flip))
        basis[:, :, 2] = broadcast(*, kernf_2f...)
        filt_onebyone!(tmp_filtered_1, image_1, kernf_2f, filter_half_size, points)
        kernf_2b = kernelfactors((gaus, gaus_der))
        filt_onebyone!(tmp_filtered_2, image_2, kernf_2b, filter_half_size, points)
        dif[:, :, 2] = tmp_filtered_2 - tmp_filtered_1

        # basis 3 - (gaus_der_flip, gaus) as forward, (gaus_der, gaus) as backward
        kernf_3f = kernelfactors((gaus_der_flip, gaus))
        basis[:, :, 3] = broadcast(*, kernf_3f...)
        filt_onebyone!(tmp_filtered_1, image_1, kernf_3f, filter_half_size, points)
        kernf_3b = kernelfactors((gaus_der, gaus))
        filt_onebyone!(tmp_filtered_2, image_2, kernf_3b, filter_half_size, points)
        dif[:, :, 3] = tmp_filtered_2 - tmp_filtered_1
    end

    dif = reshape(dif, (:, filter_num))

    # prepare matrices for linear system of equations
    A = zeros(filter_num-1, filter_num-1, pixel_count)
    b = zeros(filter_num-1, pixel_count)

    for k in 1:filter_num-1
        for l in k:filter_num-1
            @views window_sum!(A[k, l, :], dif[:, k+1] .* dif[:, l+1], image_size, window_size)
            A[l, k, :] = A[k, l, :]
        end
        @views window_sum!(b[k, :], dif[:, k+1] .* dif[:, 1] .* (-1), image_size, window_size)
    end

    # Perform Gauss elimination on all pixels in parallel:
    # coeffs will be of shape: pixel_count, filter_num-1
    @views coeffs = multi_mat_div(A, b)
    # adding ones so that all base filters have their coefficients even the first one
    all_coeffs = [ones(pixel_count) coeffs]

    # make a border mask. Is 0 if its in a border of filter_half_size size.
    window_half_size = Int64.((window_size .- 1) ./ 2)
    border_mask = parent(padarray(ones((image_size .- (2 .* window_half_size))...),
                    Fill(NaN, window_half_size, window_half_size)))
    @views all_coeffs = all_coeffs .* reshape(border_mask, (:, 1))

    k = (-filter_half_size:filter_half_size)

    # Get the displacement vector field from the filters
    u1_top = zeros(image_size);
    u1_bot = zeros(image_size);
    u2_top = zeros(image_size);
    u2_bot = zeros(image_size);

    for n in 1:filter_num
        @views u1_top[:] = u1_top[:] .- sum(transpose(basis[:, :, n]) * k) .* all_coeffs[:, n];
        @views u1_bot[:] = u1_bot[:] .+ sum(basis[:, :, n]) .* all_coeffs[:, n];

        @views u2_top[:] = u2_top[:] .- sum(basis[:, :, n] * k) .* all_coeffs[:, n];
        @views u2_bot[:] = u2_bot[:] .+ sum(basis[:, :, n]) .* all_coeffs[:, n];
    end

    u_est = 2 .* ((im .* u1_top ./ u1_bot) .+ (u2_top ./ u2_bot));

    # dont use estimations whose displacement is larger than the filter_half_size
    displacement = real(u_est).^2 + imag(u_est).^2
    displacement_mask = displacement .> filter_half_size^2
    u_est[displacement_mask] .= NaN .+ NaN .* 1im;

    return u_est, all_coeffs

end


"""
    polyfilter_lap(target::Image, source::Image; filter_num::Integer=3, max_repeats::Integer=1, display::Bool=true)

Find a transformation flow (complex displacment field), that transforms image `source` to image `target`.

# Arguments
- `target::Image`: the image we want `source` to look like.
- `source::Image`: warped image we want to transform into `target`.
- `filter_num::Integer=3`: the number of basis filters used in `single_lap` calls (so far only =3 implemented).
- `max_repeats::Integer=1`: the maximum number of times an iteration of one filter size can be repeated.
- `display::Bool=true`: use verbose prints and return an array of figures.

# Outputs
- `flow::Flow`: is the complex vector field that transforms `source` closer to `target`.
- `source_reg::Image`: is the image `source` transformed by `flow`.
- [`figs::Matrix{Figure}`: is a 2D array of `PyPlot` Figures which shows the work of the algorithm at each iteration.
For each iteration there are 3 Figures in this order: 1) current `u_est`, 2) newest addition to `u_est` `Δ_u`, 3) current `source_reg`.]

# Describtion
Implements the basic concept of `Algorithm 2` from the [paper](http://www.ee.cuhk.edu.hk/~tblu/monsite/pdfs/gilliam1701.pdf) without some features.
It uses `single_lap` iteratively; in each iteration using the transformation estimated by `single_lap` to warp the `source` image closer to
the `target` image and then using this warped closer image as the source image in the next iteration, while using progressively smaller
`filter_half_sizes` to estimate even small and faster varying displacements.

See also: [`single_lap`](@ref), [`imgshow`](@ref),[`imgshowflow`](@ref), [`warp_imgshowflow`](@ref), [`Flow`](@ref)
"""
function polyfilter_lap(target::Image, source::Image; filter_num::Integer=3, max_repeats::Integer=1, display::Bool=true)

    # convert images to floats
    source = Float32.(source)
    target = Float32.(target)

    # rescale images to have the whole [0, 1] spectrum.
    target, source = LAP_julia.rescale_intensities(target, source)
    #NOTE: a histogram match might be a good idea (https://juliaimages.org/stable/function_reference/#Images.histmatch)

    # pad with zeros if sizes difer.
    target, source = LAP_julia.pad_images(target, source)

    image_size = size(target)

    # set number of layers in the filter pyramid.
    level_num = floor(Int64, log2(minimum(size(target))/8)+1)+1

    @assert ((2^(level_num)+1) <= minimum(image_size)) "level number results in a filter larger than the size of the input images."

    # displacement init.
    u_est = zeros(size(target))

    # filter half sizes array eg. [16, 8, 4, 2, 1]
    half_size_pyramid::Array{Int64,1} = 2 .^ range(level_num-1, stop=0, length=level_num)

    # at what filter size change the interpolation strategy
    interpol_change_index = findfirst(x -> x == 2, half_size_pyramid)

    if display
        num_plots = 4
        figs = Array{Figure}(undef, level_num, max_repeats*num_plots)
    end

    source_reg = source

    for iter in 1:level_num

        if display
            println("###################")
            println("ITERATION: ", iter)
            println("filter_half_size: ", half_size_pyramid[iter])
        end

        filter_half_size = half_size_pyramid[iter]::Int
        window_size::Tuple{Int64, Int64} = 2 .* filter_half_size .* (1, 1) .+ 1
        window_half_size::Tuple{Int64, Int64} = (window_size .- 1) ./ 2

        for iter_repeat in 1:max_repeats

            Δ_u = single_lap(target, source_reg, filter_half_size, window_size, filter_num)

            # USE INPAINTING TO CORRECT U_EST:
            # 1) replicate borders
            middle_vals = Δ_u[window_half_size[1]+1:end-window_half_size[1],
                              window_half_size[2]+1:end-window_half_size[2]]
            Δ_u = parent(padarray(middle_vals, Pad(:replicate, window_half_size...)))
            # 2) inpaint middle missing values
            if all(isnan.(real(Δ_u))) # if all are NaNs
                Δ_u = zeros(image_size)
            elseif any(isnan.(Δ_u)) # if some are NaNs
                LAP_julia.inpaint.inpaint_nans!(Δ_u)
            end
            # SMOOTH U_EST WITH A GAUSSIAN FILTER:
            Δ_u = LAP_julia.clean_using_gaussain(Δ_u, window_half_size)

            # add the Δ_u to my u_est
            u_est = u_est + Δ_u

            if display
                figs[iter,1] = showflow(u_est)
                figs[iter,2] = showflow(Δ_u)
                figs[iter,3] = imgshow(source_reg)
            end

            # depending on the size of the filter used, warp the source image closer to target using u_est
            if iter < interpol_change_index
                # linear interpolation
                source_reg = LAP_julia.interpolation.warp_img(source, real(u_est), imag(u_est))
            else
                # more accurate slower interpolation TODO
                source_reg = LAP_julia.interpolation.warp_img(source, real(u_est), imag(u_est))
            end

        end

        if display
            println("###################")
        end
    end

    return u_est, source_reg, figs
end

# NOTE: can parallelize
"""
    multi_mat_div(A, b)

Return `E`, where `E[i, :]` is the solutution of the least squares problem ``\\min\\|D_ix - c_i\\|^2`` for each ``D_i``
and ``c_i``, where ``D_i`` is `A[:, :, i]`, ``c_i`` is `b[:, i]` and ``i`` is the size of the second dimension
of `b` and third dimension of `A`. In other words, each row of `E` is the solution to one matrix from `A` and it's
corresponding vector from `b`.
"""
function multi_mat_div(A, b)
    res = zeros(size(A)[1], size(A)[3])
    # Threads.@threads for j in axes(A, 4) check: (https://stackoverflow.com/questions/57678890/batch-matrix-multiplication-in-julia)
    for k in axes(A)[3]
        @views res[:, k] = qr(A[:, :, k], Val(true)) \ b[:, k]
    end
    return transpose(res)
end

# NOTE: here maybe an average could be better so the size of the window
# doest effect the number range of the results
# also might improve the speed if the pixels are saved into a tmp var
"""
    function window_sum!(filter_result, pixels, image_size, window_size)

Filter the array `pixels` of shape `image_size` with a filter that performs a sum of values of the area `window_size`
with each pixel as the center. Store the result in `filter_result`.
The filtering uses `"symmetric"` padding.
"""
function window_sum!(filter_result, pixels, image_size, window_size)
    # prepare a kernel of ones of the window size
    ones_arr_1 = ones(window_size[1])
    ones_arr_2 = ones(window_size[2])
    ones_kernel = kernelfactors((ones_arr_1, ones_arr_2))

    # filtering gets a sum of pixels of window size in each coord
    imfilter!(reshape(filter_result, image_size), reshape(pixels, image_size), ones_kernel, "symmetric")
    return nothing
end

end
