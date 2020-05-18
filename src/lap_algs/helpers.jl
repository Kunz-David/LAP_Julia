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
function multi_mat_div(A, b)
    res = zeros(size(A)[1], size(A)[3])
    # Threads.@threads for j in axes(A, 4) check: (https://stackoverflow.com/questions/57678890/batch-matrix-multiplication-in-julia)
    for k in axes(A)[3]
        @views res[:, k] = qr(A[:, :, k], Val(true)) \ b[:, k]
    end
    return transpose(res)
end


function multi_mat_div_at_points(A, b, points, image_size)
    # [2, pix_count]
    res = zeros(size(A)[1], size(A)[3])
    res = reshape(res, size(A)[1], image_size...)
    A = reshape(A, size(A)[1], size(A)[1], image_size...)
    b = reshape(b, size(b)[1], image_size...)
    for k in axes(points)[2]
        # ind = points[:, k]...
        ind = CartesianIndex(points[:, k]...)
        @views res[:, points[:, k]...] = qr(A[:, :, points[:, k]...], Val(true)) \ b[:, points[:, k]...]
    end
    res = reshape(res, size(res)[1], :)
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
    ones_arr_1 = centered(ones(window_size[1]))
    ones_arr_2 = centered(ones(window_size[2]))
    ones_kernel = kernelfactors((ones_arr_1, ones_arr_2))

    # filtering gets a sum of pixels of window size in each coord
    imfilter!(reshape(filter_result, image_size), reshape(pixels, image_size), ones_kernel, "symmetric")
    return nothing
end

function window_sum_around_points!(filter_result, pixels, image_size, window_size, points)
    # prepare a kernel of ones of the window size
    ones_arr_1 = centered(ones(window_size[1]))
    ones_arr_2 = centered(ones(window_size[2]))
    ones_kernel = kernelfactors((ones_arr_1, ones_arr_2))

    # filtering gets a sum of pixels of window size in each coord
    filt_onebyone!(reshape(filter_result, image_size...), reshape(pixels, image_size), ones_kernel, Int64((window_size[1]-1)/2), points)
    return nothing
end
