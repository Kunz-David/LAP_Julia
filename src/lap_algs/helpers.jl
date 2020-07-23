using LAP_julia
using ImageFiltering: Fill, KernelFactors.gaussian, centered, kernelfactors, imfilter!, padarray, Pad
using LinearAlgebra: qr
using PyPlot: Figure
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


# NOTE: can parallelize
"""
    multi_mat_div(A, b)

Return `E`, where `E[i, :]` is the solutution of the least squares problem ``\\min\\|D_ix - c_i\\|^2`` for each ``D_i``
and ``c_i``, where ``D_i`` is `A[:, :, i]`, ``c_i`` is `b[:, i]` and ``i`` is the size of the second dimension
of `b` and third dimension of `A`. In other words, each row of `E` is the solution to one matrix from `A` and it's
corresponding vector from `b`.
"""
function multi_mat_div_using_qr(A, b)
    res = Array{Float64,2}(undef, size(A)[1], size(A)[3])
    # Threads.@threads for j in axes(A, 4) check: (https://stackoverflow.com/questions/57678890/batch-matrix-multiplication-in-julia)
    @simd for k in axes(A)[3]
        @views res[:, k] = qr(A[:, :, k], Val(true)) \ b[:, k]
    end
    return transpose(res)
end

# TODO change to inds
function multi_mat_div_at_points(A, b, inds, image_size)
    points = inds_to_points(inds)
    # [2, pix_count]
    res = similar(A, size(A, 1) + 1, size(points, 2)) # example size = (filter_count, point_count)
    A = reshape(A, size(A)[1], size(A)[1], image_size...)
    b = reshape(b, size(b)[1], image_size...)

    for k in axes(points)[2]
        # ind = points[:, k]...
        @views res[:, points[:, k]...] .= cat(qr(A[:, :, points[:, k]...], Val(true)) \ b[:, points[:, k]...], dims=2)
    end
    res = reshape(res, size(res)[1], :)
    return transpose(res)
end

# NOTE: here maybe an average could be better so the size of the window
# doest effect the number range of the results
# also might improve the speed if the pixels are saved into a tmp var
"""
    function window_sum!(result, pixels, img_size, window)

Get a sum of values (using a "symmetric" padding on the borders) in a window `window` around each point and saving the
sum into `result`.

Note: Uses a filter algorithm with a 2D ones kernel.
Note: For small images (around `256x256`) and small filters (window of size 9 and below) is faster than [`window_sum3!`](@ref). <- factcheck this again
"""
function window_sum!(result, pixels, img_size, window)
    # prepare a kernel of ones of the window size
    ones_arr = centered(ones(window[1]))
    ones_kernel = kernelfactors((ones_arr, ones_arr))

    # filtering gets a sum of pixels of window size in each coord
    imfilter!(reshape(result, img_size), reshape(pixels, img_size), ones_kernel, "symmetric")
    return result
end

"""
    function window_sum3!(result, pixels, img_size, window)

Get a sum of values (using a "symmetric" padding on the borders) in a window `window` around each point and saving the
sum into `result`.

Note: Uses a cumsum algorithm.
Note: For small images (around `256x256`) and large filters (window of size 11 and above) is faster than [`window_sum!`](@ref). <- factcheck this again
"""
function window_sum3!(result, pixels, img_size, window)
    w = Int64.((window[1]-1)/2)
    summed = padarray(reshape(pixels, img_size), Pad(:symmetric, (w+1, w+1), (w, w)))
    @views summed[-w, :] .= 0
    @views summed[:, -w] .= 0

    summed = cumsum2d!(summed, summed)
    (c, d) = img_size

    reshape(result, img_size) .= view(summed, (1+w):(w+c), (1+w):(w+d)) .-
                                 view(summed, (1+w):(w+c), (-w):(d-w-1)) .-
                                 view(summed, (-w):(c-w-1), (1+w):(w+d)) .+
                                 view(summed, (-w):(c-w-1), (-w):(d-w-1))
    return result
end

"""
    function window_sum3!(result, pixels, img_size, window, inds)

Get a sum of values (using a "symmetric" padding on the borders) in a window `window` around `inds` and saving the
sum into `result`.

Note: Uses a cumsum algorithm.
"""
function window_sum3_at_inds!(result, pixels, img_size, window, inds)

    w = Int64.((window[1]-1)/2)
    summed = padarray(reshape(pixels, img_size), Pad(:symmetric, (w+1, w+1), (w, w)))
    @views summed[-w, :] .= 0
    @views summed[:, -w] .= 0
    summed = cumsum2d!(summed, summed)

    a = CartesianIndex(w, w)
    b = CartesianIndex(w, -w-1)
    c = CartesianIndex(-w-1, w)
    d = CartesianIndex(-w-1, -w-1)
    for (matrix_ind, result_ind)  in zip(inds, CartesianIndices(result))
        result[result_ind] = summed[matrix_ind+a] -
                          summed[matrix_ind+b] -
                          summed[matrix_ind+c] +
                          summed[matrix_ind+d]
    end
    return result
end

"""
    to_lin_index(cart_ind, size)

Get linear index from CartesianIndex `cart_ind` of 2D matrix of size `size`.
"""
@inline function to_lin_index(cart_ind, size)
    return size[2]*(cart_ind[2]-1)+cart_ind[1]
end

"""
    cumsum2d!(result, A)

Cummulative sum over 2D matrix `A`, storing the result into `result`.
"""
function cumsum2d!(result, A)
    cumsum!(result, A, dims = 1)
    cumsum!(result, result, dims = 2)
    return result
end


"""
    prepare_gaussian_filters(filter_half_size)

Get KernelFactors for the first 3 forward and backward gaussian filters from the LAP paper.
"""
function prepare_gaussian_filters(filter_half_size)

    filter_size = 2 * filter_half_size + 1

    # Calculate separable filters from basis:
    sigma = (filter_half_size + 2) / 4
    centered_inds = centered(-filter_half_size:filter_half_size)
    gaus = gaussian(sigma, filter_size)
    gaus_der = gaus .* centered_inds .* (-1)/sigma^2
    gaus_der_flip = reverse(gaus_der, dims=1)

    forward_kernels = [kernelfactors((gaus, gaus)),
                       kernelfactors((gaus, gaus_der_flip)),
                       kernelfactors((gaus_der_flip, gaus))]

    backward_kernels = [kernelfactors((gaus, gaus)),
                        kernelfactors((gaus, gaus_der)),
                        kernelfactors((gaus_der, gaus))]

    return forward_kernels, backward_kernels
end




# TODO: this can be improved by using a cumsum alg
function window_sum_around_points!(result, pixels, img_size, window, points)
    # prepare a kernel of ones of the window size
    ones_arr = centered(ones(window[1]))
    ones_kernel = kernelfactors((ones_arr, ones_arr))

    # filtering gets a sum of pixels of window size in each coord
    filt_onebyone!(reshape(result, img_size...), reshape(pixels, img_size), ones_kernel, Int64((window[1]-1)/2), points)
    # imfilter!(reshape(result, img_size), reshape(pixels, img_size), ones_kernel, "symmetric")
    return result
end

function add_at_points(A, B, inds)
    for ind in inds
        if isnan(B[ind])
            continue
        end
        @assert !isnan(A[ind])
        println(B[ind], ind)
        A[ind] = A[ind] + B[ind]
    end
    return A
end

function sum_at_points(A, inds)
    return sum([A[ind] for ind in inds])
end

"""
    gem3d!(A, b)

Reduce the linear systems of equations ``C_n x_n = d_n`` to the row echelon form. The matrixes ``C``
are slices from the first 2 dimensions of `A`, the vectors ``d`` are the slices from the first dimension
of `b` and ``n`` is the size of the third dimension of `A`.

See also: [`gem3d!`](@ref), [`back_substitution3d`](@ref)
"""
function gem3d!(A, b)
    ratios = similar(A, size(A, 3))
    for k in axes(A,1)
        for l in k+1:size(A, 2)
            ratios .= A[l, k, :] ./ A[k, k, :]
            A[l, :, :] .-= ratios' .* A[k, :, :]
            b[l, :] .-= ratios .* b[k, :]
        end
    end
end


"""
    back_substitution3d(A, b)

Solve the linear systems of equations ``C_n x_n = d_n``, where the matrixes ``C`` (in row echelon form)
are slices from the first 2 dimensions of `A`, the vectors ``d`` are the slices from the first dimension
of `b` and ``n`` is the size of the third dimension of `A`.
Returns the coefficients ``x_n`` as a 2D matrix `coeffs`, where the vectors ``x_n`` are the rows of this matrix.

See also: [`gem3d!`](@ref), [`multi_mat_div2`](@ref)
"""
function back_substitution3d(A, b)
    row_count = size(A, 1)
    coeffs = similar(A, (size(A, 3), size(A, 1)))
    for k in row_count:-1:1
        coeffs[:,k] = (b[k,:])'
        for m in (k+1):row_count
            coeffs[:,k] = coeffs[:,k]-A[k,m,:] .* coeffs[:,m]
        end
        coeffs[:,k]=coeffs[:,k]./A[k,k,:];
    end
    return coeffs
end

"""
    multi_mat_div2(A, b)

Solve the linear systems of equations ``C_n x_n = d_n``, where the matrixes ``C``
are slices from the first 2 dimensions of `A`, the vectors ``d`` are the slices from the first dimension
of `b` and ``n`` is the size of the third dimension of `A`.
Returns the coefficients ``x_n`` as a 2D matrix `coeffs`, where the vectors ``x_n`` are the rows of this matrix.

Note: Uses Gaussian elimination and back substitution.
Note: For small matrices ``C``, this is an order of magnitude faster than

See also: [`gem3d!`](@ref), [`back_substitution3d`](@ref)
"""
function multi_mat_div2(A, b)
    gem3d!(A, b)
    return back_substitution3d(A, b)
end
