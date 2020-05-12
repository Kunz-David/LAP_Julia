
##
using LAP_julia
using ComputationalResources, Images, TestImages

img = testimage("lena_gray")
img = Float32.(img)

img = gen_chess()

# genrate flow
flow = gen_rand_flow(size(img), 7, 1000);

showflow(flow)

# generate warpped image
imgw = LAP_julia.interpolation.warp_img(img, -real(flow), -imag(flow));

# params:
#lap
filter_num = 3;
filter_half_size = 10;
window_size = [21, 21];
#points
point_count = 100;
spacing = 10;

mask = padarray(ones(size(img).-(2*fhs, 2*fhs)), Fill(0, (fhs, fhs), (fhs, fhs)))
# get points
points, mag, mag_points = find_keypoints_from_gradients(img, sigma = 1, number = point_count, spacing = spacing, mask = mask);

using PyPlot
# see points
pos_x = [point.pos[1] for point in points]
pos_y = [point.pos[2] for point in points]

pos = transpose([pos_x pos_y])

points = pos

imgshow(img)
PyPlot.scatter(pos_y, pos_x, marker = :x)
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




##

# define inds
inds = (50:50, 50:50)
inds = (1:1, 50:50)

imgfilt = similar(img)

# make kernel
sigma = 3
filter_size = 21
filter_half_size = Int64((filter_size-1)/2)
gaus = ImageFiltering.KernelFactors.gaussian(sigma, filter_size)
kernel = kernelfactors((gaus, gaus))

# try the inds imfilter original
padded = padarray(img, Pad(:symmetric, filter_half_size, filter_half_size))

imgshow(padded)

imfilter!(CPU1(Algorithm.FIRTiled()), imgfilt, padded, kernel, NoPad(), inds)

for ind in CartesianIndices(img)
    ind_filt = (ind[1]:ind[1], ind[2]:ind[2])
    # println(ind_filt)
    imfilter!(CPU1(Algorithm.FIRTiled()), imgfilt, padded, kernel, NoPad(), ind_filt)
end

imgshow(imgfilt)

imgfilt2 = similar(img)

imfilter!(imgfilt2, img, kernel, "symmetric")

imgshow(imgfilt2)

imgfilt == imgfilt2


result = Array{Float32}(undef, size(points, 2))

function filt_onebyone!(imgfilt, img, kernel, fhs, points)
    padded = padarray(img, Pad(:symmetric, fhs, fhs))

    @simd for k in 1:size(points, 2)
        ind_filt = (points[1, k]-fhs:points[1, k]+fhs, points[2, k]-fhs:points[2, k]+fhs)
        imfilter!(CPU1(Algorithm.FIRTiled()), imgfilt, padded, kernel, NoPad(), ind_filt)
    end
end

imgfilt = copy(img)
imgfilt = similar(img)
imgfilt2 = similar(img)

filt_onebyone!(imgfilt, img, kernel, filter_half_size, points)
imfilter!(imgfilt2, img, kernel, "symmetric")

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








function my_imfilter!(result, image, kernel_factors, "symmetric", points)
    for k in size(points, 2)
        imfilter(points[:, k]


end
