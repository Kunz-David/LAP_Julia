
##
using LAP_julia
using ComputationalResources, Images, TestImages

img = testimage("lena_gray")
img = Float32.(img)

img = gen_chess()

# genrate flow
flow = gen_tiled_flow(size(img), 7, 1000);

showflow(flow)

# generate warpped image
imgw = LAP_julia.interpolation.warp_img(img, -real(flow), -imag(flow));

# params:
#lap
filter_num = 3;
filter_half_size = 10;
fhs = filter_half_size
window_size = [21, 21];
#points
point_count = 30;
spacing = fhs*2

mask = parent(padarray(trues(size(img).-(2*fhs, 2*fhs)), Fill(false, (fhs, fhs), (fhs, fhs))))
# get points
inds = find_edge_points(img, sigma = 1, number = point_count, spacing = spacing, mask = mask);

points = LAP_julia.inds_to_points(inds)

imgshow(img)
addpoints(inds)
gcf()


u_est, coeffs = single_lap(img, imgw, filter_half_size, window_size, 3, pos)

minimum(filter(!isnan, abs.(real(u_est))))

showflow(u_est, skip_count = 2)
showflow(u_est)

showflow(flow)

u_est = single_lap(img, imgw, filter_half_size, window_size)
LAP_julia.inpaint.inpaint_nans!(u_est)
showflow(u_est)

## profile it:


# make kernel
sigma = 3
filter_size = 21
filter_half_size = Int64((filter_size-1)/2)
#black kernel
oned = ImageFiltering.centered(ones(filter_size) .* 12123123)
kernel = kernelfactors((oned, oned))

#gaus kernel
gaus = ImageFiltering.KernelFactors.gaussian(sigma, filter_size)
kernel = kernelfactors((gaus, gaus))

# try the inds imfilter original
function filt_onebyone!(imgfilt, img, kernel, fhs, points)
    padded = padarray(img, Pad(:symmetric, fhs, fhs))

    @simd for k in 1:size(points, 2)
        ind_filt = (points[1, k]-fhs:points[1, k]+fhs, points[2, k]-fhs:points[2, k]+fhs)
        imfilter!(CPU1(Algorithm.FIRTiled()), imgfilt, padded, kernel, NoPad(), ind_filt)
    end
end

imgfilt = copy(img)
imgfilt = similar(img)
imgfilt2 = copy(img)

LAP_julia.registration.
@benchmark LAP_julia.registration.filt_onebyone!(imgfilt, img, kernel, filter_half_size, points)

@benchmark imfilter!(imgfilt2, img, kernel, "symmetric")


LAP_julia.mean(imgfilt)
LAP_julia.mean(imgfilt)

imgshow(imgfilt)
imgshow(imgfilt2)


point = [11, 11]

fhs =filter_half_size
ind_filt = (point[1]-fhs:point[1]+fhs, point[2]-fhs:point[2]+fhs)
imfilter!(CPU1(Algorithm.FIRTiled()), imgfilt, padded, kernel, NoPad(), ind_filt)
# @btime imfilter!(CPU1(Algorithm.FIRTiled()), imgfilt, padded, kernel, NoPad(), ind_filt)





time = @btime filt_onebyone!(imgfilt, img, kernel, filter_half_size, points[:, 1:20])

time = @btime imfilter!(imgfilt2, img, kernel, "symmetric")

imgfilt == imgfilt2

# compare results:
for point in points
    print(imgfilt[point] == imgfilt2[point])
    println("    ", imgfilt[point], "    ", imgfilt2[point])
end

scatter(points[1, :], points[2, :], marker=:x)
gcf()


# all_points = Array{Int64}(undef, 2, 40000)
#
# k = 1
# for ind in CartesianIndices(img)
#     all_points[:, k] = [ind[1] ind[2]]
#     global k += 1
# end

## if one by one all all_points



imgfilt = similar(img)
imgfilt2 = similar(img)

filt_onebyone!(imgfilt, img, kernel, filter_half_size, all_points)
imfilter!(imgfilt2, img, kernel, "symmetric")

@btime filt_onebyone!(imgfilt, img, kernel, filter_half_size, all_points)
@btime imfilter!(imgfilt2, img, kernel, "symmetric")

imgfilt == imgfilt2

imgshow(imgfilt)

##
using ComputationalResources


Juno.@enter imfilter!(imgfilt2, img, kernel, "symmetric")


Juno.@enter imfilter!(CPU1(Algorithm.FIRTiled()), imgfilt, padded, kernel, NoPad(), (1:1, 1:2))


##
using LAP_julia
# include("useful.jl")

img = testimage("lena_gray")
img = Float32.(img)

flow = gen_tiled_flow(size(img), 15, 100)
showflow(flow)

# generate warpped image
imgw = LAP_julia.interpolation.warp_img(img, -real(flow), -imag(flow));

# params:
#lap
filter_num = 3;
fhs = 25;
window_size = [51, 51];
#points
point_count = 25;
spacing = 40;

mask = parent(padarray(trues(size(img).-(2*fhs, 2*fhs)), Fill(false, (fhs, fhs), (fhs, fhs))))

# get points
inds = find_edge_points(img, sigma = 1, number = point_count, spacing = spacing, mask = mask);

# reform points:
pos_x = [ind[1] for ind in inds]
pos_y = [ind[2] for ind in inds]
points = LAP_julia.inds_to_points(inds)

# see points
imgshow(img)
PyPlot.scatter(pos_y, pos_x, marker = :x); gcf()

# run methods
flow_est_all = single_lap(img, imgw, fhs, window_size, 3)
# flow_est_points = single_lap(img, imgw, fhs, window_size, 3, points)
flow_est_new = single_lap_at_points(img, imgw, fhs, window_size, 3, points)

showflow(flow_est_new, disp_type=:sparse)
showflow(flow)
showflow(flow_est_new, disp_type=:sparse)

completed_flow = interpolate_flow(flow_est_new, inds)
sum(isnan.(completed_flow))

# compare resulting estimations:
# show calculated flows
showflow(completed_flow, figtitle="flow estimated from lap at points")
showflow(flow, figtitle="original flow")

# show difference and mean squared error of both
showflow(flow .- completed_flow, figtitle="orig - estim from points")
LAP_julia.mse(flow, completed_flow)
mag_mse_points = LAP_julia.vec_len(LAP_julia.mse(flow, completed_flow))
println("Points: magnitude of the vector of the mean squared error is: ", mag_mse_points)

showflow(flow .- flow_est_all, figtitle="orig - estim classic single")
inpainted = copy(flow_est_all); LAP_julia.inpaint_nans!(inpainted)
LAP_julia.mse(flow, inpainted)
mag_mse_classic = LAP_julia.vec_len(LAP_julia.mse(flow, inpainted))
println("Classic: magnitude of the vector of the mean squared error is: ", mag_mse_classic)

# whole process speed comparison:
# new method
bench_find_points = @benchmark inds = find_edge_points(img, sigma = 1, number = point_count, spacing = spacing, mask = mask)
bench_lap_points = @benchmark flow_est_new, all_coeffs = single_lap_at_points(img, imgw, fhs, window_size, 3, points)
bench_fit_points = @benchmark completed_flow = interpolate_flow(flow_est_new, inds)

# classic single lap
bench_single_lap = @benchmark flow_est_all = single_lap(img, imgw, fhs, window_size, 3)

# speedup
new_speed = median(bench_find_points.times) + median(bench_lap_points.times) + median(bench_fit_points.times)
classic_speed = median(bench_single_lap.times)
speedup = classic_speed/new_speed

println("the speedup is: " speedup)



showflow(flow_est_all)
showflow(flow_est_points)

showflow(flow)



imgshow(real(flow_est_new))
imgshow(imag(flow_est_new))

flow_points(points, flow)
flow_points(points, flow_est_new)

## add flow interpolation

meshgrid(x, y) = [repeat(x, outer=length(y)) repeat(y, inner=length(x))]
meshgrid(x::Real, y::Real) = [repeat(1:x, outer=y) repeat(1:y, inner=x)]

function test_scatter(flow_size, points, samples; method=Polyharmonic(2))
    gridPoints = meshgrid(flow_size...)'
    itp = interpolate(method, points, samples);
    interpolated = ScatteredInterpolation.evaluate(itp, gridPoints)
    gridded = reshape(interpolated, flow_size)
    return gridded
end

points[:, 1]
points[:, 1]

samples = [flow_est_new[points[:, k]...] for k in axes(points)[2]]

# multi quad
multiquad_flow_est = test_scatter(size(flow), points, samples, method=Multiquadratic(2))

showflow(multiquad_flow_est)
showflow(flow)
showflow(flow_est_new)

# thin plate
thinplate_flow_est = test_scatter(size(flow), points, samples, method=Polyharmonic(2))

showflow(thinplate_flow_est)
showflow(flow_est_new)

figure()
PyPlot.scatter([x[2] for x in inds], [x[1] for x in inds], marker = :x); gcf()


# check if it used the points unchanged
flow_points(points, flow_est_new)
flow_points(points, thinplate_flow_est) # yes
flow_points(points, multiquad_flow_est) # yes



## helper functions

flow_points(points, flow) = [flow[points[:, k]...] for k in axes(points)[2]]

showflow(flow)
showflow(flow_est_new)
showsparseflow(flow_est_new)


poly_est, source_reg, figs = polyfilter_lap(img, imgw, display=true)

showflow(poly_est, figtitle="classic poly")
showflow(flow .- poly_est, mag=1, figtitle="difference: classic poly - orig")

showflow((flow .- poly_est)[10:size(flow, 1)-10, 10:size(flow, 1)-10], mag=2)


poly_point_est, source_point_reg, figs_point = sparse_pflap(img, imgw, display=true)

poly_point_less_est, source_point_reg, figs_point, last_u = sparse_pflap(img, imgw, display=true, point_count=20, spacing=40)

showflow(poly_point_est, figtitle="points")
showflow(flow, figtitle="original")

showflow(flow .- poly_point_est, mag=1)
showflow(flow .- poly_point_less_est, mag=2)


#comapare speed:
@btime poly_point_est, source_point_reg = sparse_pflap(img, imgw, display=false)
@btime poly_est, source_reg = polyfilter_lap(img, imgw, display=false)


imgshow(imgw)
imgshow(img)

figs_point[5, 1]
figs_point[5, 2]
figs_point[7, 2]

figs_point[7, 1]
figs_point[1, 5]


## tmp hacks

function add_at_points(A, vals, inds)
    for (ind, val) in zip(inds, vals)
        A[ind] = A[ind] + val
    end
    return A
end

function add_from_matrix_at_points(A, vals, inds)
    return [val + A[ind] for (val, ind) in zip(vals, inds)]
end


A = zeros(12,12)

vals = [12, 12]
inds = [CartesianIndex(2, 3), CartesianIndex(4, 5)]

A = add_at_points(A, vals, inds)

newvals = add_from_matrix_at_points(A, vals, inds)
