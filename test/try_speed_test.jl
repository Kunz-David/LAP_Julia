### add modules
using Interpolations, ColorVectorSpace, ImageFiltering, Images
using Images
using TestImages

### setup
function imWarp(img, dx, dy)
    itp = interpolate(img, BSpline(Linear()))
    inds = indices_spatial(img)
    rng = extrema.(inds)
    imw = similar(img, eltype(itp))
    for I in CartesianIndices(inds)
        dxi, dyi = dx[I], dy[I]
        y, x = clamp(I[1]+dyi, rng[1]...), clamp(I[2]+dxi, rng[2]...)
        imw[I] = itp(y, x)
    end
    return imw
end

# get image
n = 20
mgrid = ones(512, 512);
@views mgrid[1:n:end, :] .= 0;
@views mgrid[:, 1:n:end] .= 0;
img = mgrid
# fig = PyPlot.figure(dpi = 300, figsize = (5, 5))
# PyPlot.imshow(mgrid, cmap = :gray); gcf()


#rand flow
flow = gen_tiled_flow(size(mgrid), 50);

#center flow:
ix = size(img)[1]
iy = size(img)[2]
ixhalf = Int64(ix/2)
iyhalf = Int64(iy/2)
flow = Array{ComplexF64}(undef, ix, iy)
flow[1:ixhalf, 1:iyhalf] .= -20 .+ -20 .* im;
flow[1:ixhalf, iyhalf:iy] .= 20 .+ -20 .* im;
flow[ixhalf:ix, 1:iyhalf] .= -20 .+ 20 .* im;
flow[ixhalf:ix, iyhalf:iy] .= 20 .+ 20 .* im;

dx, dy = real(flow), imag(flow);
imgw = imWarp(mgrid, dx, dy);

fig = PyPlot.figure(dpi = 300, figsize = (5, 5))
PyPlot.imshow(imgw, cmap = :gray);
ax = gca()
ax.get_ylim()
ax.set_ylim(0, size(imgw)[1])

flow[1:ixhalf, 1:iyhalf] .= 0;

showflow(flow, 20, fig=fig);

flow[1,1]
gcf()

showflow(flow, 20);
gcf()

u = gen_tiled_flow((51, 123), 10)
showflow(u)

1
LAP_julia.visualise.imgshowflow(imgw, flow)

LAP_julia.visualise.warp_imgshowflow(img, flow .* (-1))
