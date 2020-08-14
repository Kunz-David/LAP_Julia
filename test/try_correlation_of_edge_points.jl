

img, imgw = Float64.(Gray.(img)), Float64.(Gray.(imgw))

viewimg(img)
viewimg(imgw)

fhs, window = 15, [31, 31]
u_est = single_lap(img, imgw, fhs, window);

showflow(u_est)
showflow(flow)

# # show the NaN location
# nan_img = map(x -> isnan(x) ? 1.0 : 0.0, u_est);
#
# viewimg(nan_img)

function error_percentage(estim, truth)
    error = abs(estim - truth)
    return error/truth
end

function abs_error(estim, truth)
    return abs(estim - truth)
end

# good_estimation = map(x -> LAP_julia.vec_len(error_percentage(x...)), zip(u_est, flow.*(-1)))
# maximum(filter(!isnan, good_estimation))
# minimum(filter(!isnan, good_estimation))
#
# good_estimation_r = maximum(filter(!isnan, good_estimation)) .- good_estimation
#
#
# maximum(filter(!isnan, good_estimation_r))
# minimum(filter(!isnan, good_estimation_r))
#
# figure()
#
# imgshow(img)
# imshow(good_estimation, cmap = :RdYlGn_r, alpha = 0.3); gcf()
# colorbar(ax=gca()); gcf()


nccs=zeros(100)

for k in 1:100
    img, imgw, flow = gen_init(:lena, :uniform, flow_args=[rand(Complex{Float64}), 5]);
    grad, imge = LAP_julia.gradient_magnitude(img);
    flow_est = single_lap(img, imgw, 10, [21,21]);
    percentage_errors = map(x -> LAP_julia.vec_len(error_percentage(x...)), zip(flow_est, flow.*(-1)));
    nccs[k] = ncc(imge[.!isnan.(flow_est)], percentage_errors[.!isnan.(flow_est)])
end

figure()
plt.hist(nccs, 20); gcf()

mean(nccs)


## with mri data
nccs=zeros(100)

for k in 1:100
    img, imgw, flow = gen_init(:lena, :uniform, flow_args=[rand(Complex{Float64}), 5]);
    grad, imge = LAP_julia.gradient_magnitude(img);
    flow_est = single_lap(img, imgw, 10, [21,21]);
    percentage_errors = map(x -> LAP_julia.vec_len(error_percentage(x...)), zip(flow_est, flow.*(-1)));
    nccs[k] = ncc(imge[.!isnan.(flow_est)], percentage_errors[.!isnan.(flow_est)])
end

figure()
plt.hist(nccs, 20); gcf()
base_path = "/Users/MrTrololord/Downloads/mri_head_homography/"
target_images_path = joinpath(base_path, "images", "target")

function gen_mri(path_folder, flow_function)
    img = load(joinpath(path_folder, "img$(string(rand((collect(1:60))))).png"))
    flow = flow_function(size(img))
    imgw = warp_img(img, -real.(flow), -imag.(flow))
    return img, imgw, flow
end

flow_function(img_size) = gen_uniform_flow(img_size, rand(Complex{Float64}), 5)

img, imgw, flow = gen_mri(target_images_path, flow_function)

# test on mri
nccs=zeros(100)

for k in 1:100
    img, imgw, flow = gen_mri(target_images_path, flow_function);
    grad, imge = LAP_julia.gradient_magnitude(imgw);
    img, imgw = Float64.(Gray.(img)), Float64.(Gray.(imgw))
    flow_est = single_lap(img, imgw, 10, [21,21]);
    percentage_errors = map(x -> LAP_julia.vec_len(error_percentage(x...)), zip(flow_est, flow.*(-1)));
    nccs[k] = ncc(imge[.!isnan.(flow_est)], percentage_errors[.!isnan.(flow_est)])
end

fig, ax = subplots(dpi = 400)
plt.hist(nccs, 20);
title("NCC of Displacement Error and Magnitude of Gradient");
savefig("../plots/ncc_displacement_error_and_grad.png")

gcf()
mean(nccs)

# abs

abs_err = map(x -> LAP_julia.vec_len(abs_error(x...)), zip(u_est, flow.*(-1)))
maximum(filter(!isnan, abs_err))
minimum(filter(!isnan, abs_err))

imgshow(img)
imshow(abs_err, cmap = :RdYlGn_r, alpha = 0.3); gcf()
colorbar(ax=gca()); gcf()




vec, edge = LAP_julia.gradient_magnitude(img)

maximum(filter(!isnan, 1 ./good_estimation))
minimum(filter(!isnan, 1 ./good_estimation))

filter(!isnan, normalize_to_zero_one(1 ./good_estimation))

minimum(filter(!isnan, normalize_to_zero_one(1 ./good_estimation)))

using LAP_julia: normalize_to_zero_one

imgshow(map(x-> isnan(x) ? 0 : x, normalize_to_zero_one(1 ./good_estimation)))

collect(colorview(RGB, normalize_to_zero_one(edge), normalize_to_zero_one(edge), normalize_to_zero_one(1 ./good_estimation)))

collect(colorview(Gray, img.*(1/extrema(img)[2]))) .+ collect(colorview(RGBA, good_estimation, zeros(size(good_estimation)), zeros(size(good_estimation)), fill(0.3, size(good_estimation)...)))


imgoverlay_v2(img, good_estimation)

imgshow(ones(112,112))
