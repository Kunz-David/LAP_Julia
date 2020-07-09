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

#FIXME: use window_size not fhs
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
    function window_sum!(result, pixels, img_size, window)

Filter the array `pixels` of shape `img_size` with a filter that performs a sum of values of the area `window`
with each pixel as the center. Store the result in `result`.
The filtering uses `"symmetric"` padding.

Note: For small images (around `256x256`) and small filters (window of size 9 and below) is faster than [`window_sum3!`](@ref).
"""
function window_sum!(result, pixels, img_size, window)
    # prepare a kernel of ones of the window size
    ones_arr = centered(ones(window))
    ones_kernel = kernelfactors((ones_arr, ones_arr))

    # filtering gets a sum of pixels of window size in each coord
    imfilter!(reshape(result, img_size), reshape(pixels, img_size), ones_kernel, "symmetric")
    return nothing
end

"""
    function window_sum3!(result, pixels, img_size, window)

Filter the array `pixels` of shape `img_size` with a filter that performs a sum of values of the area `window`
with each pixel as the center. Store the result in `result`.
The filtering uses `"symmetric"` padding.

Note: For small images (around `256x256`) and large filters (window of size 11 and above) is faster than [`window_sum!`](@ref).
"""
function window_sum3!(result, pixels, img_size, w)
    summed = padarray(reshape(pixels, img_size), Pad(:symmetric, (w, w)))
    summed = cumsum2d!(summed, summed)
    B = padarray(summed, Fill(0, (1, 1), (0, 0)))
    c = img_size[1]

    reshape(result, img_size) .= view(B, (1+w):(w+c), (1+w):(w+c)) .-
              view(B, (1+w):(w+c), (-w):(c-w-1)) .-
              view(B, (-w):(c-w-1), (1+w):(w+c)) .+
              view(B, (-w):(c-w-1), (-w):(c-w-1))
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


# TODO: this can be improved by using a cumsum alg
function window_sum_around_points!(filter_result, pixels, image_size, window_size, points)
    # prepare a kernel of ones of the window size
    ones_arr_1 = centered(ones(window_size[1]))
    ones_arr_2 = centered(ones(window_size[2]))
    ones_kernel = kernelfactors((ones_arr_1, ones_arr_2))

    # filtering gets a sum of pixels of window size in each coord
    filt_onebyone!(reshape(filter_result, image_size...), reshape(pixels, image_size), ones_kernel, Int64((window_size[1]-1)/2), points)
    return nothing
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
    ratios = zeros(size(A, 3))
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
    coeffs = zeros(size(A, 3), size(A, 1))
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
