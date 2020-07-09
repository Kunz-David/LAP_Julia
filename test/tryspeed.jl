include("../src/lap_algs/helpers.jl")

k, l = 1, 1

A = rand(123, 123, 123)

var1 = A[l, :, :]
time1 = @benchmark A[l, :, :] = rand(123,123)
time1 = mean(time1.times)


B = rand(12345, 4, 4)

var2 = A[:, l, :]
time2 = @benchmark A[:, l, :] = rand(123,123)
time2 = mean(time2.times)


var3 = A[:, :, l]
time3 = @benchmark A[:, :, l] = rand(123,123)
time3 = mean(time3.times)


##
A = rand(2, 2, 65536);
@benchmark A[k, l, :] .= rand(65536)

A = rand(65536, 2, 2);
@benchmark A[:, k, l] .= rand(65536)

##


N = 5
C = rand(N, N, 12345)
d = rand(N, 12345)

time4 = @benchmark multi_mat_div2(C, d)
time4 = mean(time4.times)

time5 = @benchmark multi_mat_div(C, d)
time5 = mean(time5.times)

time5/time4

out4 = multi_mat_div2(C, d)
out5 = multi_mat_div(C, d)

maximum(abs.(out4-out5)) < 10e-7


## test single lap with and without gem

img, imgw, flow = gen_init();
fhs = 23
window = [47, 47]

showflow(flow)

# no gem
oldflow = LAP_julia.classic_alg(img, imgw, fhs, window)
time6 = @benchmark LAP_julia.classic_alg(img, imgw, fhs, window)
time6 = mean(time6.times)

# gem
newflow = LAP_julia.classic_alg(img, imgw, fhs, window)
time7 = @benchmark LAP_julia.classic_alg(img, imgw, fhs, window)
time7 = mean(time7.times)

time6/time7


## multi mat div 2    with views and without

@benchmark (@views multi_mat_div2($C, $d))
@benchmark multi_mat_div2($C, $d)

#--> views are better here

## try cumsum of matrix
img = rand(2565, 2565)

function cumsum2d(A)
    result = cumsum(A, dims = 1)
    cumsum!(result, result, dims = 2)
    return result
end

function cumsum2d!(result, A)
    cumsum!(result, A, dims = 1)
    cumsum!(result, result, dims = 2)
    return result
end

window_sum2
S(x,y) = a(x,y) + S(x-1,y) + S(x,y-1) - S(x-1,y-1)

include("../src/lap_algs/helpers.jl")

imgg = copy(img)
window_sum!(imgg, img, (256, 256), (47,47))

imgg

##

small = [1 2 3; 4 5 6; 7 8 9]

@inline function get_sum_around_ind(ind, B, radius) # square window
    hw = radius
    x, y = ind[1], ind[2]
    return B[x + hw, y + hw] - B[x + hw, y - hw - 1] - B[x - hw - 1, y + hw] + B[x - hw - 1, y - hw - 1]
end

function window_sum2!(result, pixels, img_size, window, timer)
    @timeit timer "borders" begin
        w = Int64((window-1)/2)
        summed = padarray(pixels, Pad(:symmetric, (w, w),
                                           (w, w)))
    end
    @timeit timer "cumsum" begin
        summed = cumsum2d!(summed, summed)
        B = padarray(summed, Fill(0, (1, 1), (1, 1)))
    end

    @timeit timer "for loop" begin
        @simd for ind in CartesianIndices(pixels)
            result[ind] = get_sum_around_ind(ind, B, w)
        end
    end
    return result
end


function window_sum3!(result, pixels, img_size, window, timer)
    @timeit timer "borders" begin
        w = Int64((window-1)/2)
        summed = padarray(pixels, Pad(:symmetric, (w, w),
                                           (w, w)))
    end
    @timeit timer "cumsum" begin
        summed = cumsum2d!(summed, summed)
        B = padarray(summed, Fill(0, (1, 1), (0, 0)))
    end
    c = size(pixels, 1)

    @timeit timer "calculations" begin
        result .= view(B, (1+w):(w+c), (1+w):(w+c)) .-
                  view(B, (1+w):(w+c), (-w):(c-w-1)) .-
                  view(B, (-w):(c-w-1), (1+w):(w+c)) .+
                  view(B, (-w):(c-w-1), (-w):(c-w-1))
    end
end

A = rand(234, 234)
result3 = Float64.(copy(A))
result2 = Float64.(copy(A))
result1 = Float64.(copy(A))


# timed
timer = TimerOutput();

reset_timer!(timer)
@timeit timer "whole" window_sum3!(result3, A, size(A), 21, timer)
print_timer(timer)

reset_timer!(timer)
@timeit timer "whole" window_sum2!(result2, A, size(A), 21, timer)
print_timer(timer)


result3
result2

A = small
# bench
@benchmark window_sum3!(result3, A, size(A), 9, timer)
@benchmark window_sum3!(result3, A, size(A), 9)
@benchmark window_sum2!(result2, A, size(A), 31, timer)
@benchmark window_sum!(result1, A, size(A), (9,9))

maximum(abs.(result1 - result2))

## time classic alg improved
TimerOutputs.enable_debug_timings(LAP_julia.lap)

img, imgw, flow = gen_init()

fhs = 23
window = [47, 47]

fhs = 3
window = [7, 7]
newflow = LAP_julia.classic_alg(img, imgw, fhs, window)
