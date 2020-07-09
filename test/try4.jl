### chessboard

tile_size = 50
board_size = 4 # must be even
mini_board = [zeros(tile_size, tile_size) ones(tile_size, tile_size);
              ones(tile_size, tile_size) zeros(tile_size, tile_size)]

chessboard = repeat(mini_board, outer=(convert(Integer,(board_size/2)), convert(Integer,(board_size/2))))
img = chessboard

img = gen_chess();

###

f(x :: Array{<:Real, 2}) = println(x)

Array{Float64, 2} <: Array{Real, 2}

Float64 <: Real

function g()
    ar::Array{Real, 2} = ones(1,1)
end

g()

f(ones(1, 1))


###

using TestImages
using PyPlot

img = testimage("mandril_gray");
img = reverse(img, dims=1)
img = chessboard

point_count = 100;
spacing = 30;

coverage = 1;

(pi * spacing^2 * point_count)/length(img)

r = sqrt(( * ))

points, mag, mag_points = LAP_julia.gradient_points.find_keypoints_from_gradients(img, sigma = 5, number = point_count, spacing = spacing)

imgshow(mag_points)
imgshow(mag)

pos_x = [point.pos[1] for point in points]
pos_y = [point.pos[2] for point in points]

imgshow(img)
PyPlot.scatter(pos_y, pos_x, marker = :x)
gcf()

### probability
using PyPlot

img = testimage("mandril_gray");
img = chessboard

# all(vec(reshape(img, (1, :))) == img[:])

points, mag, mag_points = LAP_julia.gradient_points.find_keypoints_from_gradients_p_field(img, sigma = 1, number = 100, spacing = 30)

wsample(Array(1:2), [-0, 1], 12)

imgshow(mag_points)
imgshow(mag)

pos_x = [point.pos[1] for point in points]
pos_y = [point.pos[2] for point in points]

imgshow(img)
PyPlot.scatter(pos_y, pos_x, marker = :x)
gcf()



###
using LAP_julia

img = gen_chess()

# genrate flow
flow = gen_tiled_flow(size(img), 20, 1000);

showflow(flow)

# generate warpped image
imgw = LAP_julia.interpolation.warp_img(img, real(flow), imag(flow));

# params:
#lap
image_1 = img;
image_2 = imgw;
filter_num = 3;
filter_half_size = 32;
window_size = [65, 65];
#points
point_count = 100;
spacing = 10;

# get points
points, mag, mag_points = find_keypoints_from_gradients(img, sigma = 1, number = point_count, spacing = spacing);

using PyPlot
# see points
pos_x = [point.pos[1] for point in points]
pos_y = [point.pos[2] for point in points]

pos = transpose([pos_x pos_y])
pos[:, 1]

imgshow(img)
PyPlot.scatter(pos_y, pos_x, marker = :x)
gcf()


u_est, coeffs = single_lap(image_1, image_2, filter_half_size, window_size, 3, pos);

showflow(u_est)


### try imfilter with indices:
img = chessboard;

filter_size = 21
filter_half_size = (filter_size - 1) / 2
# Calculate separable filters from basis:
sigma = (filter_half_size + 2) / 4
gaus = ImageFiltering.KernelFactors.gaussian(sigma, filter_size)

using PyPlot
using Plots
# 1
kernel = kernelfactors((gaus, gaus))
kernp1 = broadcast(*, kernel...)
surface([-filter_half_size:filter_half_size], -filter_half_size:filter_half_size, kernp1)

@time imgfilt1 = ImageFiltering.imfilter(chessboard, kernel, "symmetric")
imgshow(imgfilt1)
surface(imgfilt1)

# @time imgfilt1 = ImageFiltering.imfilter(chessboard, kernel, "symmetric", inds=[1, 1])

imgfilt = similar(img)

inds = CartesianIndex(50, 50)

imfilter!(CPU1(Algorithm.FIRTiled()), imgfilt, img, kernel, "symmetric", inds)

imfilter!(CPU1(Algorithm.FIR()), imgfilt, img, (), NoPad(), inds)

imgshow(imgfilt)
surface(imgfilt1)


### exapmle
using ImageFiltering, Images, ComputationalResources

# Create a sample image:
tile_size = 50
board_size = 4 # must be even
mini_board = [zeros(tile_size, tile_size) ones(tile_size, tile_size);
              ones(tile_size, tile_size) zeros(tile_size, tile_size)]
chessboard = repeat(mini_board, outer=(convert(Integer,(board_size/2)), convert(Integer,(board_size/2))))
img = chessboard

# Create a filter:
sigma = 3
filter_size = 21
gaus = ImageFiltering.KernelFactors.gaussian(sigma, filter_size)
kernel = kernelfactors((gaus, gaus))



typeof(kernel) <: ImageFiltering.ProcessedKernel
typeof(kernel) <: Union{ImageFiltering.ArrayLike, ImageFiltering.Laplacian}

# params of imfilter function
imgfilt = copy(img)
# inds = CartesianIndex(50, 50)

# The filtering that errors:
# imfilter!(CPU1(Algorithm.FIRTiled()), imgfilt, img, kernel, "symmetric", (30:40, 30:40))
imfilter!(CPU1(Algorithm.FIRTiled()), imgfilt, img, kernel, Pad(:Symmetric), (50:50, 50:50))

imgshow(imgfilt)


### alternative to imfilter inds
function f(img_pix)
    sum(img_pix .* [0 1/10 0;
                    1/10 6/10 1/10;
                    0 1/10 0])
end

mapped_img = mapwindow(f, img, (3, 3), "symmetric", (2:8, 2:8))

imgshow(mapped_img)

###

flow = gen_tiled_flow(size(img), 30, 1000)



one = [1 + 0im, 2 - im]
two = [0 + 6im, 3 + im]

one = [1 + 0im]
two = [0 + im]

import LAP_julia: helpers.angle_rms

angle_rms(one, two)

rad2deg.(angle.(one).-angle.(two))


### mse


function mse2(x, y)
    @assert (size(x) == size(y)) "sizes of $x and $y don't match"
    return 1/length(x) * sum((x .- y).^2)
end

function mse1(x, y)
    @assert (size(x) == size(y)) "sizes of $x and $y don't match"
    return mean((x .- y).^2)
end

function mean(x)
    sum(x)/length(x)
end

ran = rand(123)
dan = rand(123)

mse1(ran, dan)
mse2(ran, dan)
mse1(ran, dan) ≈ mse2(ran, dan)


####
