### TRY numero 2

tile_size = 4
board_size = 8 # must be even
mini_board = [zeros(tile_size, tile_size) ones(tile_size, tile_size);
              ones(tile_size, tile_size) zeros(tile_size, tile_size)]

chessboard = repeat(mini_board, outer=(convert(Integer,(board_size/2)), convert(Integer,(board_size/2))))
chessboard = repeat(mini_board, outer=(2, 2))

img = chessboard



### test clean for u


chess_norm = chessboard
chess_norm_f = Float32.(chess_norm)

chess_warp_f = copy(chess_norm_f)
for k in 2:size(chess_warp_f)[1]
    for l in 2:size(chess_warp_f)[2]
        chess_warp_f[k-1, l-1] = chess_warp_f[k, l]
        if l == size(chess_warp_f)[2]
            chess_warp_f[k-1, l] = chess_warp_f[k, l]
        end
        if k == size(chess_warp_f)[1]
            chess_warp_f[k, l-1] = chess_warp_f[k, l]
        end
    end
end

filter_half_size = 2
filter_size = filter_half_size * 2 + 1

u_est, coeffs = LAP_julia.lap.single_lap(chess_norm_f, chess_warp_f, 3, filter_half_size, [filter_size, filter_size])

window_half_size = [1, 1]
u_out = LAP_julia.helpers.smooth_with_gaussian(u_est[3:end-2, 3:end-2], window_half_size)

LAP_julia.visualise.showflow(u_est[3:end-2, 3:end-2])

LAP_julia.visualise.showflow(u_out)

### test gen rand flow
for k in 1:2
    out = gen_tiled_flow((100,100),10)
    showflow(out, 1)
end

out = gen_tiled_flow((100,100),10)
showflow(out)
# --->works

### Tim Holy INTERPOLATION example

using Interpolations, ColorVectorSpace, ImageFiltering, Images

function imWarp(img, dx, dy)
    itp = interpolate(img, BSpline(Linear()))
    inds = indices_spatial(img)
    rng = extrema.(inds)
    imw = similar(img, eltype(itp))
    println(rng[1])
    for I in CartesianIndices(inds)
        dxi, dyi = dx[I], dy[I]
        y, x = I[1]+dyi, I[2]+dxi
        if is_in(y, rng[1]) && is_in(x,rng[2])
            imw[I] = itp(y, x)
        else
            imw[I] = 0
        end
    end
    return imw
end

function is_in(x, (low, high))
    return x >= low && x <= high
end

is_in(1.0, (1, 512))

# Usage demo
using TestImages
img = testimage("mandril");
img

# Let's use a random warp. We blur it with a Gaussian kernel to make it reasonably smooth.
using ImageFiltering

#test1:
kern = KernelFactors.IIRGaussian((10,10))  # IIRGaussian is fast even for very large σ
dx, dy = imfilter(100*randn(size(img)), kern), imfilter(100*randn(size(img)), kern);
imgw = imWarp(img, dx, dy);
imgw

#test2:
flow = gen_tiled_flow(size(img), 50);
dx, dy = real(flow), imag(flow);
imgw = imWarp(img, dx, dy);
imgw

#test3:
n = 20
grid = ones(512, 512);
@views grid[1:n:end, :] .= 0;
@views grid[:, 1:n:end] .= 0;
fig = PyPlot.figure(dpi = 400, figsize = (5, 5))
PyPlot.imshow(grid, cmap = :gray); gcf()
flow = gen_tiled_flow(size(grid), 50);
dx, dy = real(flow), imag(flow);
imgw = imWarp(grid, dx, dy);

fig = PyPlot.figure(dpi = 400, figsize = (5, 5))
PyPlot.imshow(imgw, cmap = :gray); gcf()


fig = PyPlot.figure(dpi = 400, figsize = (5, 5))
PyPlot.close_figs()
showflow(flow,200)

## pyplot disp image !!!! Does not work for TestImages

PyPlot.imshow(chessboard, cmap = :gray)
gcf()


###
using Interpolations

t = 0:.1:1
x = sin.(2π*t)
y = cos.(2π*t)
A = hcat(x,y)

itp = Interpolations.scale(interpolate(A, (BSpline(Cubic(Natural(OnCell()))), NoInterp())), t, 1:2)

tfine = 0:.01:1
xs, ys = [itp(t,1) for t in tfine], [itp(t,2) for t in tfine];


using PyPlot

figure()
scatter(x, y, label="knots")
PyPlot.plot(xs, ys, label="spline")
gcf()

###
using Interpolations
using PyPlot

t = 0:.1:1
x = sin.(2π*t)

itp = Interpolations.interpolate(x, (BSpline(Quadratic(Natural(OnGrid())))))
itp = Interpolations.interpolate(x, (BSpline(Linear())));
sitp = scale(itp, t);
#sitp = itp

tfine = 0:.01:1;
xs = [sitp(t) for t in tfine]


figure()
scatter(t, x, label="knots")
PyPlot.plot(tfine, xs, label="spline")
gcf()

### get rid of unwanted image replication on the right and bottom of images in the interpolation

img = testimage("mandril_gray");
ix = size(img)[1]
iy = size(img)[2]

ixhalf = Int64(ix/2)
iyhalf = Int64(iy/2)

# flow all to centere
flow = Array{ComplexF64}(undef, ix, iy)
flow[1:ixhalf, 1:iyhalf] .= -20 .+ -20 .* im;
flow[1:ixhalf, iyhalf:iy] .= 20 .+ -20 .* im;
flow[ixhalf:ix, 1:iyhalf] .= -20 .+ 20 .* im;
flow[ixhalf:ix, iyhalf:iy] .= 20 .+ 20 .* im;

flow = LAP_julia.helpers.smooth_with_gaussian(flow, 60);

showflow(flow)

dx, dy = real(flow), imag(flow);
imgw = imWarp(img, dx, dy)
imgw


save("img.png", imgw)

###
