## speedup window sum at inds

image_size = size(img)
pixels = length(img)

result = zeros(pixels)
pixels = reshape(img, pixels)
window = [25, 25]
inds = rand(CartesianIndices(img), 50)


function window_sum3_at_inds!(result, pixels, img_size, window, inds)
    w = Int64.((window[1]-1)/2)
    summed = padarray(reshape(pixels, img_size), Pad(:symmetric, (w, w)))
    summed = cumsum2d!(summed, summed)
    B = padarray(summed, Fill(0, (1, 1), (0, 0)))
    c = img_size[1]

    reshape(result, img_size) .= view(B, (1+w):(w+c), (1+w):(w+c)) .-
                                 view(B, (1+w):(w+c), (-w):(c-w-1)) .-
                                 view(B, (-w):(c-w-1), (1+w):(w+c)) .+
                                 view(B, (-w):(c-w-1), (-w):(c-w-1))
end


## whole image
import LAP_julia: window_sum3!, cumsum2d!, window_sum!

img_size = size(img)
pixel_count = length(img)

result = zeros(pixel_count)
result1 = zeros(pixel_count)
pixels = reshape(img, pixels)
window = [25, 25]
inds = rand(CartesianIndices(img), 50)

@benchmark window_sum3!(result, pixels, img_size, window, timer)
@benchmark window_sum3_new_pad!(result1, pixels, img_size, window, timer)

result == result1

timer = TimerOutput("timer")
@timeit timer "whole" begin
    B = window_sum3_new_pad!(result, pixels, img_size, window, timer)
end
print_timer(timer)


## attempt at speedup

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
        res[result_ind] = summed[matrix_ind+a] -
                          summed[matrix_ind+b] -
                          summed[matrix_ind+c] +
                          summed[matrix_ind+d]
    end
    return res
end

res_inds = deepcopy(result)
res = deepcopy(result)

res = window_sum3!(res, pixels, img_size, window)
res_inds = window_sum3_at_inds!(res_inds, pixels, img_size, window, CartesianIndices(img_size))


@benchmark window_sum3!(res, pixels, img_size, window)
@benchmark window_sum3_at_inds!(res_inds, pixels, img_size, window, rand(CartesianIndices(img_size), 300))


result2 = deepcopy(result1)


res = window_sum3!(result1, pixels, img_size, window)
res_old = window_sum!(result2, pixels, img_size, window)

res_new_pad == res


B1 == B

result == result1

A = reshape(1:90000, 300, 300)

out = padarray(A, Pad(:symmetric, (2, 2), (1, 1)))



@benchmark out[:, 1] .= 0
@benchmark @views out[:, 1] .= 0

@benchmark out[1, :] .= 0
@benchmark @views out[1, :] .= 0



t = randn(5,5,2)
[ sum( t[ 1,1,: ] ), sum( t[ 5,5,: ] ) ]

sum.()




sum((x->x[I]).([A,B,C,D]))
