
### TEST pflap

## using grid

# generate image GRID
n = 20
mgrid = ones(512, 512);
@views mgrid[1:n:end, :] .= 0;
@views mgrid[:, 1:n:end] .= 0;
img = mgrid

# genrate flow
flow = gen_tiled_flow(size(img), 50);

# generate warpped image
imgw = LAP_julia.interpolation.imWarp(img, real(flow), imag(flow))


## using mandrill image

# genrate image MANDRILL
img = testimage("mandril_gray")
imgr = reverse(img, dims=1)
img = imgr

# genrate flow
flow = gen_tiled_flow(size(img), 30);
showflow(flow)

# generate warpped image
imgw = LAP_julia.interpolation.imWarp(img, real(flow), imag(flow));

# show imgw
imgshow(Float32.(imgw))

imgshow(img)

## using chessboard image

tile_size = 5
board_size = 4 # must be even
mini_board = [zeros(tile_size, tile_size) ones(tile_size, tile_size);
              ones(tile_size, tile_size) zeros(tile_size, tile_size)]
chessboard = repeat(mini_board, outer=(convert(Integer,(board_size/2)), convert(Integer,(board_size/2))))
img = chessboard

img = gen_chess()

# genrate flow
flow = gen_tiled_flow(size(img), 20, 1000);

# generate warpped image
imgw = LAP_julia.interpolation.warp_img(img, real(flow), imag(flow))

# show imgw
imgshow(imgw)

# show flow
showflow(flow)

maximum(filter(!isnan, real(flow)))
minimum(filter(!isnan, real(flow)))

### RUN OPTIFLOW:

u_est, source_reg, figs = pflap(img, imgw);

# compare with a single lap:
u_sin_est = single_lap(img, imgw, 32, [65, 65], 3)
showflow(u_sin_est)

LAP_julia.inpaint.inpaint_nans!(u_sin_est)
single_img = LAP_julia.interpolation.imWarp_replicate(imgw, real(u_sin_est), imag(u_sin_est))
imgshow(single_img)
imgshow(imgw)

### compare results

# estimation
maximum(abs.(imag(u_est)))
showflow(u_est)
showflow(u_est .* (-1) .- flow)

# rewarp check
dewarp = LAP_julia.interpolation.imWarp(imgw, real(u_est), imag(u_est))
imgshow(Float32.(dewarp))
imgshow(Float32.(img))

imgshow(source_reg)

# truth
maximum(abs.(imag(flow)))
sum(abs.(real(flow)))/length(real(flow))
showflow(flow, mag = 1)

# show source_reg
fig = PyPlot.figure(dpi = 300, figsize = (5, 5))
PyPlot.imshow(source_reg, cmap = :gray);
gcf()

# show imgw
fig = PyPlot.figure(dpi = 300, figsize = (5, 5))
PyPlot.imshow(Float64.(imgw), cmap = :gray);
gcf()

















# k = (-1:1) * ones(1,3)
# A =      [1     2     3;
#           4     5     6;
#           7     8     9];
#
# transpose(A[:]) * k[:]
#
# sum(transpose(A) * (-1:1))
#
# # l0
#
# l = transpose((-1:1) * ones(1,3))
#
# transpose(A[:]) * l[:]
#
# sum(A * (-1:1))
1
