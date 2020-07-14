
using TimerOutputs
using ImageFiltering: Fill, KernelFactors.gaussian, centered, kernelfactors, imfilter!, padarray, Pad
using LinearAlgebra: qr
using PyPlot: Figure

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

function goo()
    to = TimerOutput()
    @timeit_debug to "sleepytime" begin
        println(getfield(lap, :timeit_debug_enabled)())
        println("hello")
        sleep(1)
    end
    print_timer(to)
    return nothing
end

"""
    single_lap(image_1::Image, image_2::Image, filter_half_size::Integer, window, filter_num::Integer=3)

Return an estimate of a smoothly varying flow of size of `image_1` which is the displacement that transforms `image_2` closer to `image_1`.

# Arguments:
- `image_1::Image`: grayscale image 1.
- `image_2::Image`: grayscale image 2.
- `filter_half_size`: the half size of the base of the gaussian filters used.
- `window`: the size of the local window (list of 2 ints) usually same as filter_size.
- `filter_num::Integer=3`: the number of basis filters used (so far only =3 implemented).

See also: [`polyfilter_lap`](@ref), [`showflow`](@ref), [`warp_imgshowflow`](@ref), [`imgshowflow`](@ref)
"""
function single_lap(image_1::Image,
                    image_2::Image,
                    filter_half_size::Integer,
                    window,
                    filter_num::Integer=3;
                    timer::TimerOutput=TimerOutput("LAP"),
                    display::Bool=false)


    pixel_count = length(image_1)
    image_size = size(image_1)
    filter_size = 2 * filter_half_size + 1

    # prepare 2D filter basis for later use:
    basis = similar(image_1, (filter_size, filter_size, filter_num))

    # dif = (image forward - image backward) of each filter
    dif = similar(image_1, (image_size..., filter_num))

    @timeit_debug timer "filtering" begin
        if filter_num == 3
            # get the relevant forward and backward filters
            forward_ker, backward_ker = prepare_gaussian_filters(filter_half_size)

            # temporary place to store filtered images
            tmp_filtered_1 = similar(image_1)
            tmp_filtered_2 = similar(image_2)

            # basis 1 - (gaus, gaus) both forward and backward
            basis[:, :, 1] = broadcast(*, forward_ker[1]...)
            imfilter!(tmp_filtered_1, image_1, forward_ker[1], "symmetric")
            imfilter!(tmp_filtered_2, image_2, backward_ker[1], "symmetric")
            dif[:, :, 1] = tmp_filtered_2 - tmp_filtered_1

            # basis 2 - (gaus, gaus_der_flip) as forward, (gaus, gaus_der) as backward
            basis[:, :, 2] = broadcast(*, forward_ker[2]...)
            imfilter!(tmp_filtered_1, image_1, forward_ker[2], "symmetric")
            imfilter!(tmp_filtered_2, image_2, backward_ker[2], "symmetric")
            dif[:, :, 2] = tmp_filtered_2 - tmp_filtered_1

            # basis 3 - (gaus_der_flip, gaus) as forward, (gaus_der, gaus) as backward
            basis[:, :, 3] = broadcast(*, forward_ker[3]...)
            imfilter!(tmp_filtered_1, image_1, forward_ker[3], "symmetric")
            imfilter!(tmp_filtered_2, image_2, backward_ker[3], "symmetric")
            dif[:, :, 3] = tmp_filtered_2 - tmp_filtered_1
        end
    end

    dif = reshape(dif, (:, filter_num))

    # prepare matrices for linear system of equations
    A = similar(image_1, (filter_num-1, filter_num-1, pixel_count))
    b = similar(image_1, (filter_num-1, pixel_count))

    @timeit_debug timer "prepare A and b" begin
        w = Int64.((window[1]-1)/2)
        for k in 1:filter_num-1
            for l in k:filter_num-1
                @timeit_debug timer "window sum part 1" @views window_sum3!(A[k, l, :], dif[:, k+1] .* dif[:, l+1], image_size, window)
                A[l, k, :] = A[k, l, :]
            end
            @timeit_debug timer "window sum part 2" @views window_sum3!(b[k, :], dif[:, k+1] .* dif[:, 1] .* (-1), image_size, window)
        end
    end

    # Perform Gauss elimination on all pixels in parallel:
    # coeffs will be of shape: pixel_count, filter_num-1
    @timeit_debug timer "multli mat div 2" begin
        @views coeffs = multi_mat_div2(A, b)
    end
    # adding ones so that all base filters have their coefficients even the first one
    all_coeffs = [ones(pixel_count) coeffs]

    # make a border mask. Is 0 if its in a border of filter_half_size size.
    window_half_size = Int64.((window .- 1) ./ 2)
    border_mask = parent(padarray(ones((image_size .- (2 .* window_half_size))...),
                    Fill(NaN, window_half_size, window_half_size)))
    @views all_coeffs = all_coeffs .* reshape(border_mask, (:, 1))

    k = (-filter_half_size:filter_half_size)

    # Get the displacement vector field from the filters
    u1_top = similar(image_1);
    u1_bot = similar(image_1);
    u2_top = similar(image_1);
    u2_bot = similar(image_1);

    @timeit_debug timer "calculate flow" begin
        for n in 1:filter_num
            @views u1_top[:] = u1_top[:] .- sum(transpose(basis[:, :, n]) * k) .* all_coeffs[:, n];
            @views u1_bot[:] = u1_bot[:] .+ sum(basis[:, :, n]) .* all_coeffs[:, n];

            @views u2_top[:] = u2_top[:] .- sum(basis[:, :, n] * k) .* all_coeffs[:, n];
            @views u2_bot[:] = u2_bot[:] .+ sum(basis[:, :, n]) .* all_coeffs[:, n];
        end

        u_est = 2 .* ((im .* u1_top ./ u1_bot) .+ (u2_top ./ u2_bot));
    end

    # dont use estimations whose displacement is larger than the filter_half_size
    displacement = real(u_est).^2 .+ imag(u_est).^2
    displacement_mask = displacement .> filter_half_size^2
    u_est[displacement_mask] .= NaN .+ NaN .* 1im;

    if (display && getfield(LAP_julia, :timeit_debug_enabled)())
        print_timer(timer)
        println()
    end
    return u_est
end



function single_lap_at_points(image_1::Image,
                              image_2::Image,
                              filter_half_size::Integer,
                              window,
                              points,
                              filter_num::Integer=3;
                              timer::TimerOutput=TimerOutput("Sparse LAP"))

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
    basis = similar(image_1, (filter_size, filter_size, filter_num))

    # dif = (points of forward image - points of backward image) for each filter
    dif = similar(image_1, (image_size..., filter_num))

    @timeit_debug timer "filtering" begin
        if filter_num == 3
            # temporary place to store filtered images
            tmp_filtered_1 = fill(NaN, image_size)
            tmp_filtered_2 = fill(NaN, image_size)

            # basis 1 - (gaus, gaus) both forward and backward
            kernf_1 = kernelfactors((gaus, gaus))
            basis[:, :, 1] = broadcast(*, kernf_1...)
            # filt_onebyone!(tmp_filtered_1, image_1, kernf_1, filter_half_size, points)
            # filt_onebyone!(tmp_filtered_2, image_2, kernf_1, filter_half_size, points)
            imfilter!(tmp_filtered_1, image_1, kernf_1, "symmetric")
            imfilter!(tmp_filtered_2, image_2, kernf_1, "symmetric")
            dif[:, :, 1] = tmp_filtered_2 - tmp_filtered_1

            # basis 2 - (gaus, gaus_der_flip) as forward, (gaus, gaus_der) as backward
            kernf_2f = kernelfactors((gaus, gaus_der_flip))
            basis[:, :, 2] = broadcast(*, kernf_2f...)
            # filt_onebyone!(tmp_filtered_1, image_1, kernf_2f, filter_half_size, points)
            imfilter!(tmp_filtered_1, image_1, kernf_2f, "symmetric")
            kernf_2b = kernelfactors((gaus, gaus_der))
            # filt_onebyone!(tmp_filtered_2, image_2, kernf_2b, filter_half_size, points)
            imfilter!(tmp_filtered_2, image_2, kernf_2b, "symmetric")
            dif[:, :, 2] = tmp_filtered_2 - tmp_filtered_1

            # basis 3 - (gaus_der_flip, gaus) as forward, (gaus_der, gaus) as backward
            kernf_3f = kernelfactors((gaus_der_flip, gaus))
            basis[:, :, 3] = broadcast(*, kernf_3f...)
            # filt_onebyone!(tmp_filtered_1, image_1, kernf_3f, filter_half_size, points)
            imfilter!(tmp_filtered_1, image_1, kernf_3f, "symmetric")
            kernf_3b = kernelfactors((gaus_der, gaus))
            # filt_onebyone!(tmp_filtered_2, image_2, kernf_3b, filter_half_size, points)
            imfilter!(tmp_filtered_2, image_2, kernf_3b, "symmetric")
            dif[:, :, 3] = tmp_filtered_2 - tmp_filtered_1
        end
    end # "filtering"

    dif = reshape(dif, (:, filter_num))

    # prepare matrices for linear system of equations
    # A = zeros(filter_num-1, filter_num-1, pixel_count)
    # b = zeros(filter_num-1, pixel_count)
    A = fill(NaN, filter_num-1, filter_num-1, pixel_count)
    b = fill(NaN, filter_num-1, pixel_count)

    @timeit_debug timer "prepare A and b" begin
        # TODO: check whether window_sum_around_points calculates only around points and not on whole window.
        for k in 1:filter_num-1
            for l in k:filter_num-1
                @timeit_debug timer "window sum part 1" begin
                    # @views window_sum_around_points!(A[k, l, :], dif[:, k+1] .* dif[:, l+1], image_size, window, points)
                    @views window_sum3!(A[k, l, :], dif[:, k+1] .* dif[:, l+1], image_size, window)
                    A[l, k, :] = A[k, l, :]
                end
            end
            @timeit_debug timer "window sum part 2" begin
                # @views window_sum_around_points!(b[k, :], dif[:, k+1] .* dif[:, 1] .* (-1), image_size, window, points)
                @views window_sum3!(b[k, :], dif[:, k+1] .* dif[:, 1] .* (-1), image_size, window)
            end
        end
    end # "prepare A and b"

    # Perform Gauss elimination on all pixels in parallel:
    # coeffs will be of shape: pixel_count, filter_num-1
    @timeit_debug timer "multi mat div" begin
        @views coeffs = multi_mat_div_at_points(A, b, points, image_size)
    end
    # adding ones so that all base filters have their coefficients even the first one
    all_coeffs = [ones(pixel_count) coeffs]
    println("zeros in coeffs: ", count(x -> x==0, coeffs), "non zeros: ", count(x -> x !=0, coeffs))

    # make a border mask. Is 0 if its in a border of filter_half_size size.
    window_half_size = Int64.((window .- 1) ./ 2)
    border_mask = parent(padarray(ones((image_size .- (2 .* window_half_size))...),
                    Fill(NaN, window_half_size, window_half_size)))
    @views all_coeffs = all_coeffs .* reshape(border_mask, (:, 1))

    k = (-filter_half_size:filter_half_size)

    # Get the displacement vector field from the filters
    u1_top = zeros(image_size);
    u1_bot = zeros(image_size);
    u2_top = zeros(image_size);
    u2_bot = zeros(image_size);

    @timeit_debug timer "calculate flow" begin
        for n in 1:filter_num
            @views u1_top[:] .= u1_top[:] .- sum(transpose(basis[:, :, n]) * k) .* all_coeffs[:, n];
            @views u1_bot[:] .= u1_bot[:] .+ sum(basis[:, :, n]) .* all_coeffs[:, n];

            @views u2_top[:] .= u2_top[:] .- sum(basis[:, :, n] * k) .* all_coeffs[:, n];
            @views u2_bot[:] .= u2_bot[:] .+ sum(basis[:, :, n]) .* all_coeffs[:, n];
        end

        u_est = 2 .* ((im .* u1_top ./ u1_bot) .+ (u2_top ./ u2_bot));
    end # "calculate flow"

    # dont use estimations whose displacement is larger than the filter_half_size
    displacement = real(u_est).^2 + imag(u_est).^2
    displacement_mask = displacement .> filter_half_size^2
    u_est[displacement_mask] .= NaN .+ NaN .* 1im;

    return u_est

end


function single_lap_at_points_win_sum1(image_1::Image,
                              image_2::Image,
                              filter_half_size::Integer,
                              window,
                              points,
                              filter_num::Integer=3;
                              timer::TimerOutput=TimerOutput("Sparse LAP"))

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
    basis = similar(image_1, (filter_size, filter_size, filter_num))

    # dif = (points of forward image - points of backward image) for each filter
    dif = similar(image_1, (image_size..., filter_num))

    @timeit_debug timer "filtering" begin
        if filter_num == 3
            # temporary place to store filtered images
            tmp_filtered_1 = fill(NaN, image_size)
            tmp_filtered_2 = fill(NaN, image_size)

            # basis 1 - (gaus, gaus) both forward and backward
            kernf_1 = kernelfactors((gaus, gaus))
            basis[:, :, 1] = broadcast(*, kernf_1...)
            # filt_onebyone!(tmp_filtered_1, image_1, kernf_1, filter_half_size, points)
            # filt_onebyone!(tmp_filtered_2, image_2, kernf_1, filter_half_size, points)
            imfilter!(tmp_filtered_1, image_1, kernf_1, "symmetric")
            imfilter!(tmp_filtered_2, image_2, kernf_1, "symmetric")
            dif[:, :, 1] = tmp_filtered_2 - tmp_filtered_1

            # basis 2 - (gaus, gaus_der_flip) as forward, (gaus, gaus_der) as backward
            kernf_2f = kernelfactors((gaus, gaus_der_flip))
            basis[:, :, 2] = broadcast(*, kernf_2f...)
            # filt_onebyone!(tmp_filtered_1, image_1, kernf_2f, filter_half_size, points)
            imfilter!(tmp_filtered_1, image_1, kernf_2f, "symmetric")
            kernf_2b = kernelfactors((gaus, gaus_der))
            # filt_onebyone!(tmp_filtered_2, image_2, kernf_2b, filter_half_size, points)
            imfilter!(tmp_filtered_2, image_2, kernf_2b, "symmetric")
            dif[:, :, 2] = tmp_filtered_2 - tmp_filtered_1

            # basis 3 - (gaus_der_flip, gaus) as forward, (gaus_der, gaus) as backward
            kernf_3f = kernelfactors((gaus_der_flip, gaus))
            basis[:, :, 3] = broadcast(*, kernf_3f...)
            # filt_onebyone!(tmp_filtered_1, image_1, kernf_3f, filter_half_size, points)
            imfilter!(tmp_filtered_1, image_1, kernf_3f, "symmetric")
            kernf_3b = kernelfactors((gaus_der, gaus))
            # filt_onebyone!(tmp_filtered_2, image_2, kernf_3b, filter_half_size, points)
            imfilter!(tmp_filtered_2, image_2, kernf_3b, "symmetric")
            dif[:, :, 3] = tmp_filtered_2 - tmp_filtered_1
        end
    end # "filtering"

    dif = reshape(dif, (:, filter_num))

    # prepare matrices for linear system of equations
    # A = zeros(filter_num-1, filter_num-1, pixel_count)
    # b = zeros(filter_num-1, pixel_count)
    A = fill(NaN, filter_num-1, filter_num-1, pixel_count)
    b = fill(NaN, filter_num-1, pixel_count)

    @timeit_debug timer "prepare A and b" begin
        # TODO: check whether window_sum_around_points calculates only around points and not on whole window.
        for k in 1:filter_num-1
            for l in k:filter_num-1
                @timeit_debug timer "window sum part 1" begin
                    # @views window_sum_around_points!(A[k, l, :], dif[:, k+1] .* dif[:, l+1], image_size, window, points)
                    @views window_sum!(A[k, l, :], dif[:, k+1] .* dif[:, l+1], image_size, window)
                    A[l, k, :] = A[k, l, :]
                end
            end
            @timeit_debug timer "window sum part 2" begin
                # @views window_sum_around_points!(b[k, :], dif[:, k+1] .* dif[:, 1] .* (-1), image_size, window, points)
                @views window_sum!(b[k, :], dif[:, k+1] .* dif[:, 1] .* (-1), image_size, window)
            end
        end
    end # "prepare A and b"

    # Perform Gauss elimination on all pixels in parallel:
    # coeffs will be of shape: pixel_count, filter_num-1
    @timeit_debug timer "multi mat div" begin
        @views coeffs = multi_mat_div_at_points(A, b, points, image_size)
    end
    # adding ones so that all base filters have their coefficients even the first one
    all_coeffs = [ones(pixel_count) coeffs]

    # make a border mask. Is 0 if its in a border of filter_half_size size.
    window_half_size = Int64.((window .- 1) ./ 2)
    border_mask = parent(padarray(ones((image_size .- (2 .* window_half_size))...),
                    Fill(NaN, window_half_size, window_half_size)))
    @views all_coeffs = all_coeffs .* reshape(border_mask, (:, 1))

    k = (-filter_half_size:filter_half_size)

    # Get the displacement vector field from the filters
    u1_top = zeros(image_size);
    u1_bot = zeros(image_size);
    u2_top = zeros(image_size);
    u2_bot = zeros(image_size);

    @timeit_debug timer "calculate flow" begin
        for n in 1:filter_num
            @views u1_top[:] = u1_top[:] .- sum(transpose(basis[:, :, n]) * k) .* all_coeffs[:, n];
            @views u1_bot[:] = u1_bot[:] .+ sum(basis[:, :, n]) .* all_coeffs[:, n];

            @views u2_top[:] = u2_top[:] .- sum(basis[:, :, n] * k) .* all_coeffs[:, n];
            @views u2_bot[:] = u2_bot[:] .+ sum(basis[:, :, n]) .* all_coeffs[:, n];
        end

        u_est = 2 .* ((im .* u1_top ./ u1_bot) .+ (u2_top ./ u2_bot));
    end # "calculate flow"

    # dont use estimations whose displacement is larger than the filter_half_size
    displacement = real(u_est).^2 + imag(u_est).^2
    displacement_mask = displacement .> filter_half_size^2
    u_est[displacement_mask] .= NaN .+ NaN .* 1im;

    return u_est

end
