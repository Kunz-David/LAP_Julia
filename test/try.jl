
#--- load modules
using .LAP_julia
using Test, Images, Colors
using FileIO
using CSV
using Plots, ImageView: imshow

#--- load random gray ANHIR image
base_path = "/Users/MrTrololord/Google_Drive/cvut/bakalarka/anhir/"
dataset_path = "/Users/MrTrololord/Google_Drive/cvut/bakalarka/anhir/dataset/"
landmark_path = "/Users/MrTrololord/Google_Drive/cvut/bakalarka/anhir/landmarks/"
loc_table = CSV.read(base_path * "location_table.csv")

train_rows = loc_table[loc_table[:status] .== "training", :]

const TEST_SIZE = 1

random_rows = train_rows[rand(1:size(train_rows, 1), TEST_SIZE), :]

row_used = (random_rows[1, :])[Symbol("Column1")]

source_loc = (random_rows[1, :])[Symbol("Source image")]
target_loc = (random_rows[1, :])[Symbol("Target image")]
source_in = load(dataset_path * source_loc)
target_in = load(dataset_path * target_loc)

gray_source = Gray.(source_in)
gray_target = Gray.(target_in)

#--- load chessboard image1

tile_size = 5
board_size = 8 # must be even
mini_board = [zeros(tile_size, tile_size) ones(tile_size, tile_size);
              ones(tile_size, tile_size) zeros(tile_size, tile_size)]

chessboard = repeat(mini_board, outer=(convert(Integer,(board_size/2)), convert(Integer,(board_size/2))))
chessboard = repeat(mini_board, outer=(2, 2))

imshow(chessboard)

# imshow(chessboard);

#---
small_gray_source = imresize(gray_source, ratio = 0.1)

# LAP filters:
filter_size = 21
filter_half_size = (filter_size - 1) / 2
use_one_dim_filters = false
# Calculate separable filters from basis:
sigma = (filter_half_size + 2) / 4
centered_inds = ImageFiltering.centered(-filter_half_size:filter_half_size)
gaus = ImageFiltering.KernelFactors.gaussian((sigma, sigma), (filter_size, filter_size))
gaus = ImageFiltering.KernelFactors.gaussian(sigma, filter_size)
gaus_der = gaus .* centered_inds .* (-1)/sigma^2
gaus_der_flip = reverse(gaus_der, dims=1)

# kernel basis:
# 1
kernf1 = kernelfactors((gaus, gaus))
kernp1 = broadcast(*, kernf1...)
surface(-filter_half_size:filter_half_size, -filter_half_size:filter_half_size, kernp1)

# 2
# forward
kernf2 = kernelfactors((gaus, gaus_der_flip))
kernp2 = broadcast(*, kernf2...)
surface(-filter_half_size:filter_half_size, -filter_half_size:filter_half_size, kernp2)
# backward
kernf2 = kernelfactors((gaus, gaus_der))
kernp2 = broadcast(*, kernf2...)
surface(-filter_half_size:filter_half_size, -filter_half_size:filter_half_size, kernp2)

# 3
# forward
kernf3 = kernelfactors((gaus_der_flip, gaus))
kernp3 = broadcast(*, kernf3...)
surface(-filter_half_size:filter_half_size, -filter_half_size:filter_half_size, kernp3)
# backward
kernf3 = kernelfactors((gaus_der, gaus))
kernp3 = broadcast(*, kernf3...)
surface(-filter_half_size:filter_half_size, -filter_half_size:filter_half_size, kernp3)




surface(-filter_half_size:filter_half_size, -filter_half_size:filter_half_size, kernp)


#--- application to images

# kernels:
# 1
@time imgfilt1 = ImageFiltering.imfilter(chessboard, kernf1, "symmetric")
imshow(imgfilt1)
surface(imgfilt1)

# 2
imgfilt2 = ImageFiltering.imfilter(chessboard, kernf2, "symmetric")
imshow(imgfilt2)
surface(imgfilt2)

# 3
imgfilt3 = ImageFiltering.imfilter(chessboard, kernf3, "symmetric")
imshow(imgfilt3)
surface(imgfilt3)

#--- imfilter with prealloc

imgfilt1_prealloc = similar(chessboard)

@time imgfilt1_tmp = ImageFiltering.imfilter!(imgfilt1_prealloc, chessboard, kernf1, "symmetric")
imshow(imgfilt1_prealloc)
surface(imgfilt1_prealloc)

#---


dif = zeros(4, 4, 3)

dif[:, :, 2] = 1:16
dif[:, :, 3] = ones(4, 4).+1

dif

dif = reshape(dif, (:, 3))

dif

#---

A = zeros(2, 2, 3)
b = zeros(2, 3)

A = rand(100, 100, 100)
b = rand(100, 100)

A[:, :, 1] = [1 0; 0 1]
A[:, :, 2] = [1 0; 0 1]
A[:, :, 3] = [1 0; 0 1]

b[:, 1] = [1; 1]
b[:, 2] = [0; 2]
b[:, 3] = [0; 4]

A[:, :, 1]\b[:, 1]

# ??? try it with * broadcast
res = (\).(A, b)

tmp = (*).(A, b)

save = zeros(100, 100)

@time multi_mat_div!(save, A, b)

save

function multi_mat_div!(save, A, b)
    for k in axes(A)[3]
        @views save[:, k] = A[:, :, k] \ b[:, k]
    end
end


function multi_mat_div(A, b)
    res = zeros(size(A)[1], size(A)[3])
    # Threads.@threads for j in axes(A, 4) check: (https://stackoverflow.com/questions/57678890/batch-matrix-multiplication-in-julia)
    for k in axes(A)[3]
        @views res[:, k] = A[:, :, k] \ b[:, k]
    end
    return res
end

@time multi_mat_div(A, b)

#---

one = ones(12);
two = reshape(one, (3, 4))

one[:] = two

#--- pad with zeros

using ImageFiltering

image_size = size(chessboard)
filter_half_size = 5;

ones(image_size.-(2*filter_half_size))

border_mask = parent(padarray(ones(image_size.-(2*filter_half_size)), Fill(NaN, (filter_half_size, filter_half_size),(filter_half_size, filter_half_size))))

new_board = chessboard .* border_mask

imshow(new_board)

all_coeffs = ones(broadcast((*), image_size...), 3)

new_coeffs = all_coeffs .* reshape(border_mask, (:, 1))

imshow(reshape(new_coeffs[:, 1], image_size))

imshow(reshape(chessboard[border .== 1.0] .= NaN, (170, 170)))


#---

using BenchmarkTools
using LinearAlgebra

A = rand(11, 11);

A[1, :] = ones(11)

k = (-filter_half_size:filter_half_size)

k1 = repeat(k, outer=11)

sum(transpose(A) * k)

dot(A[:], repeat(k, outer=11))

B = A[:]


#---




u_est, coeffs = single_lap(chessboard, chessboard, 3, 15, [15, 15])


#---

ones_arr_1 = ones(window_size[1])
ones_arr_2 = ones(window_size[2])
ones_kernel = kernelfactors((ones_arr_1, ones_arr_2))
# filtering gets a sum of pixels of window size in each coord

pixels = reshape(chessboard, (:, 1))
image_size = size(chessboard)
window_size = [15, 15]

filter_result = similar(chessboard)

ImageFiltering.imfilter!(filter_result, reshape(pixels, image_size), shit, "symmetric")

imshow(filter_result)

#---

include("../src/lap.jl")

lap.window_sum!(filter_result, pixels, image_size, window_size)

imshow(filter_result)

#---

uv = zeros(size(chessboard)..., 2)
uv[:, :, 1] = real(u_est)
uv[:, :, 2] = imag(u_est)

Plots.quiver([1,2,3], [1,2,3], quiver=([NaN,1,1],[NaN,2,3]))



u_est, coeffs = single_lap(chessboard, chessboard, 3, 15, [15, 15])

p1 = figure()
uv_flow = zeros(size(u_est)..., 2)
uv_flow[:, :, 1] = real(u_est)
uv_flow[:, :, 2] = imag(u_est)

skip_count = 20
n = skip_count

trimmed_real = zeros(size(uv_flow[:,:,1]))
trimmed_real[1:n:end, 1:n:end] = uv_flow[1:n:end, 1:n:end, 1]
trimmed_imag = zeros(size(uv_flow[:,:,2]))
trimmed_imag[1:n:end, 1:n:end] = uv_flow[1:n:end, 1:n:end, 2]

PyPlot.quiver(trimmed_real, trimmed_imag)
gcf()


#---
using FileIO: load

board_path = "/Users/MrTrololord/Google_Drive/cvut/bakalarka/PolyFilter_LAP/test_boards/"

chess_norm = load(board_path * "chessboard.png")
chess_warp = load(board_path * "chess_warp.png")


chess_norm_f = Float32.(chess_norm)
chess_warp_f = Float32.(chess_warp)



u_est, coeffs = single_lap(chess_norm_f, chess_warp_f, 3, 15, [15, 15])

using LinearAlgebra: qr

qr([1 1; 1 1], Val(true)) \ [1; 1]

#--- TEST single_lap
using PyPlot: quiver, gcf, figure

# *************************
# 1) MATLAB RANDOM FLOW
# *************************
board_path = "/Users/MrTrololord/Google_Drive/cvut/bakalarka/PolyFilter_LAP/test_boards/"

# load chessboards from matlab generated flow
chess_norm = load(board_path * "chessboard.png")
chess_warp = load(board_path * "chess_warp.png")

chess_norm_f = Float32.(chess_norm)
chess_warp_f = Float32.(chess_warp)


u_est, coeffs = LAP_julia.lap.single_lap(chess_norm_f, chess_warp_f, 3, 15, [15, 15])

# *************************
# 2) ONE PIXEL IMAGE SHIFT
# *************************

chess_norm = chessboard
chess_norm_f = Float32.(chess_norm)

chess_warp_f = copy(chess_norm_f)
for k in 2:size(chess_warp_f)[1]
    for l in 2:size(chess_warp_f)[2]
        chess_warp_f[k-1, l-1] = chess_warp_f[k, l]
    end
end

filter_half_size = 2
filter_size = filter_half_size * 2 + 1

u_est, coeffs = LAP_julia.lap.single_lap(chess_norm_f, chess_warp_f, 3, filter_size, [filter_size, filter_size])

u_est[10, 10]
real(u_est)
imag(u_est)

p1 = figure()
quiver(real(u_est), imag(u_est))
gcf()
