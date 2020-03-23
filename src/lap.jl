module lap
export single_lap, polyfilter_lap

using ImageFiltering: Fill, KernelFactors.gaussian, centered, kernelfactors, imfilter!, padarray
using LinearAlgebra: qr

function test()
    print("testing r   eviseakdkd")
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
    single_lap(image_1, image_2, base_filters, )

# input:
- `image_1` ... grayscale image 1
- `image_2` ... grayscale image 2
- `filter_num` ... number of basis filters used (so far only =3 implemented)
- `filter_size` ... filter has size: `filter_size` \times `filter_size`
- `window_size` ... size of local window (list of 2 ints)


"""
function single_lap(image_1, image_2, filter_num, filter_size, window_size)

    print("testing attention please.")
    image_size = size(image_1)
    filter_half_size = round(Int, (filter_size - 1) / 2)



    # Prepare filters:
    # 1D
    use_one_dim_filters = true

    # Calculate separable filters from basis:
    sigma = (filter_half_size + 2) / 4
    centered_inds = centered(-filter_half_size:filter_half_size)
    gaus = gaussian(sigma, filter_size)
    gaus_der = gaus .* centered_inds .* (-1)/sigma^2
    gaus_der_flip = reverse(gaus_der, dims=1)

    # OLD way:
    # gaussian(x,s) = exp.(-x.*x/2/s^2)
    # K0=ceil([K,L]); --> this is filter half size
    # Gx = gaus(-filter_half_size:filter_half_size, sigma)
    # Gdx = (-filter_half_size:filter_half_size) .* Gx
    # Gx = Gx ./ sum(Gx)
    # Gdx = Gdx ./ sum((-filter_half_size:filter_half_size) .* Gdx)
    # Gdix = reverse(Gdx, dims=1)

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
    A = zeros(filter_num-1, filter_num-1, prod(image_size))
    b = zeros(filter_num-1, prod(image_size))

    for k in 1:filter_num-1
        for l in k:filter_num-1
            @views window_sum!(A[k, l, :], dif[:, k+1] .* dif[:, l+1], image_size, window_size)
            A[l, k, :] = A[k, l, :]
        end
        @views window_sum!(b[k, :], dif[:, k+1] .* dif[:, 1] .* (-1), image_size, window_size)
    end

    # Perform Gauss elimination on all pixels in parallel:
    # coeffs will be of shape: prod(image_size), filter_num-1
    @views coeffs = multi_mat_div(A, b)
    # adding ones so that all base filters have their coefficients even the first one
    all_coeffs = [ones(prod(image_size)) coeffs]

    # make a border mask. Is 0 if its in a border of filter_half_size size.
    border_mask = parent(padarray(ones(image_size.-(2*filter_half_size)), Fill(NaN, (filter_half_size, filter_half_size),(filter_half_size, filter_half_size))))
    @views all_coeffs = all_coeffs .* reshape(border_mask, (:, 1))

    k = (-filter_half_size:filter_half_size)

    # Get the displacement vector field from the filters
    u1_top = zeros(image_size);
    u1_bot = zeros(image_size);
    u2_top = zeros(image_size);
    u2_bot = zeros(image_size);

    for n in 1:filter_num
        @views u1_top[:] = u1_top[:] - sum(transpose(basis[:, :, n]) * k) .* all_coeffs[:, n];
        @views u1_bot[:] = u1_bot[:] + sum(basis[:, :, n]) .* all_coeffs[:, n];

        @views u2_top[:] = u2_top[:] - sum(basis[:, :, n] * k) .* all_coeffs[:, n];
        @views u2_bot[:] = u2_bot[:] + sum(basis[:, :, n]) .* all_coeffs[:, n];
    end

    u_est = 2 .* ((u1_top ./ u1_bot) .+ (im .* u2_top ./ u2_bot));

    # TODO: debug this part
    # dont use estimations whose displacement is larger than the filter_half_size
    # displacement_mask = real(u_est).^2 + imag(u_est).^2
    # u_est[displacement_mask .> filter_half_size^2] .= NaN;

    return u_est, all_coeffs

end


function polyfilter_lap(target, source)

    # choose filter basis size. 3 or 6
    filter_num = 3

    # rescale images to have the whole [0, 1] spectrum.
    target, source = rescale_intensities(target, source)

    # pad with zeros if sizes difer.
    target, source = pad_images(target, source)

    image_size = size(target)

    # set number of layers in the filter pyramid.
    level_num = floor(log2(minimum(size(target))/8)+1)

    # displacement init.
    u_est = zeros(size(target))

    @assert ((2^(level_num+1)+1) <= minimum(image_size)) "level number results in a filter larger than the size of the input images."

end

# NOTE: can parallelize
function multi_mat_div(A, b)
    res = zeros(size(A)[1], size(A)[3])
    # Threads.@threads for j in axes(A, 4) check: (https://stackoverflow.com/questions/57678890/batch-matrix-multiplication-in-julia)
    for k in axes(A)[3]
        @views res[:, k] = qr(A[:, :, k], Val(true)) \ b[:, k]
        qr([1 1; 1 1], Val(true))
    end
    return transpose(res)
end

# NOTE: here maybe an average could be better so the size of the window
# doest effect the number range of the results
# also might improve the speed if the pixels are saved into a tmp var
function window_sum!(filter_result, pixels, image_size, window_size)
    # prepare a kernel of ones of the window size
    ones_arr_1 = ones(window_size[1])
    ones_arr_2 = ones(window_size[2])
    ones_kernel = kernelfactors((ones_arr_1, ones_arr_2))

    # filtering gets a sum of pixels of window size in each coord
    imfilter!(reshape(filter_result, image_size), reshape(pixels, image_size), ones_kernel, "symmetric")
    #return filtered_result[:]
end

end
