
viewimg(img) = colorview(Gray, img)

plots_folder = "/Users/MrTrololord/Google_Drive/cvut/bakalarka/plots"
base_path = "/Users/MrTrololord/Downloads/head_mri/"
target_images_path = joinpath(base_path, "images", "target")
source_images_path = joinpath(base_path, "images", "source")
k = 22

img = load(joinpath(target_images_path, "img$(k).png"))
img = img.*(1/extrema(img)[2])

imgw = load(joinpath(source_images_path, "img$(k).png"))
imgw = imgw.*(1/extrema(imgw)[2])

save(joinpath(plots_folder, "head_mri_target.png"), img)
save(joinpath(plots_folder, "head_mri_source.png"), imgw)

over = imgoverlay_v2(img, imgw)
save(joinpath(plots_folder, "head_mri_overlay.png"), over)


# flow = gen_uniform_flow(size(img), rand(ComplexF64), 10)
# flow = gen_quad_flow(size(img), 30)
imgw = warp_img(img, -real.(flow), -imag.(flow))
# imgww = LAP_julia.smooth_with_gaussian!(imgw, 1)

imgshow(imgw)

COUNT = 500
max_disp = 40
flow_gen_func(img_size) = gen_homo_flow(img_size, max_disp)
cur_flow = flow_gen_func(size(img))
source = warp_img(img, -real.(cur_flow), -imag.(cur_flow))


save(joinpath(plots_folder, "head_mri_target.png"), img)
save(joinpath(plots_folder, "head_mri_source_other.png"), source)

img
source


showflow(imresize(cur_flow, ratio=0.1))

flip_out = reverse(cur_flow, dims = 1)
flip_out = real.(flip_out) .- imag(flip_out) .* 1im

showflow(imresize(flip_out.*(-1), ratio=0.1))


method = sparse_pflap
timer = TimerOutput(string(method));
method_kwargs = Dict(:display => false, :timer => timer, :match_source_histogram => true, :point_count => 600)
flow_est, source_reg, timer, time_in_secs = test_registration_alg(method, img, source, cur_flow, timer=method_kwargs[:timer], method_kwargs=method_kwargs);

flip_est = reverse(flow_est, dims = 1)
flip_est = real.(flip_est) .- imag(flip_est) .* 1im


showflow(flip_est, figtitle="Estimated Deformation")
head_mri_flow_est_path = joinpath(plots_folder, "head_mri_flow_est.png")
savefig(head_mri_flow_est_path)


showflow(flip_out.*(-1), figtitle="Ground Truth Deformation")
head_mri_flow_path = joinpath(plots_folder, "head_mri_flow.png")
savefig(head_mri_flow_path)






img, imgw = Float64.(Gray.(img)), Float64.(Gray.(imgw))

fhs, window = 15, [31, 31]
u_est = single_lap(img, imgw, fhs, window);

showflow(u_est)
showflow(flow)


# show the NaN location
nan_img = map(x -> isnan(x) ? 1.0 : 0.0, u_est);

viewimg(nan_img)

function error_percentage(estim, truth)
    error = abs(estim - truth)
    return error/truth
end

function abs_error(estim, truth)
    return abs(estim - truth)
end

good_estimation = map(x -> LAP_julia.vec_len(error_percentage(x...)), zip(u_est, flow.*(-1)))
maximum(filter(!isnan, good_estimation))
minimum(filter(!isnan, good_estimation))

good_estimation_r = maximum(filter(!isnan, good_estimation)) .- good_estimation


maximum(filter(!isnan, good_estimation_r))
minimum(filter(!isnan, good_estimation_r))

figure()

imgshow(img)
imshow(good_estimation, cmap = :RdYlGn_r, alpha = 0.3); gcf()
colorbar(ax=gca()); gcf()


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








method = sparse_pflap_psnr
timer = TimerOutput(string(method));
method_kwargs = Dict(:display => false, :timer => timer, :match_source_histogram => true)
flow_est, source_reg, timer, time_in_secs = time_reg_alg(method, img, imgw, timer=method_kwargs[:timer], method_kwargs=method_kwargs);


source_reg_zeros = warp_img(imgw, -real.(flow_est), -imag.(flow_est), border_strat=:zeros)
img


## see points

imgshow(imgw, figtitle = "Source")
imgshow(source_reg, figtitle = "Registered Source")
imgshow(img, figtitle = "Target")


imgoverlay_v2(img, source_reg_zeros)
imgoverlay_v2(img, imgw)



showflow(flow.*(-1), figtitle = "Truth Flow")
showflow(flow_est, figtitle = "Estimated Flow")


inds, mag = find_edge_points(img, debug=true)

imgshow(mag)
addpoints(inds)


method = sparse_pflap_psnr
timer = TimerOutput(string(method));
method_kwargs = Dict(:display => false, :timer => timer, :match_source_histogram => true, :point_count => 600)
flow_est, source_reg, timer, time_in_secs = time_reg_alg(method, img, imgw, timer=method_kwargs[:timer], method_kwargs=method_kwargs);




method = pflap
timer = TimerOutput(string(method));
method_kwargs = Dict(:display => false, :timer => timer, :match_source_histogram => true, :prefilter => false)
flow_est, source_reg, timer, time_in_secs = time_reg_alg(method, img, imgw, timer=method_kwargs[:timer], method_kwargs=method_kwargs);



------------------------------------------------------------
                                              Time
                                      ----------------------
           Tot / % measured:               4.34s / 100%

Section                       ncalls     time   %tot     avg
------------------------------------------------------------
sparse_pflap_psnr                  1    4.33s   100%   4.33s
  single filter pyramid level      8    4.16s  95.9%   519ms
    sparse lap                    12    3.85s  88.7%   320ms
      filtering                   12    3.47s  80.2%   289ms
      prepare A and b             12    370ms  8.54%  30.8ms
        window sum part 1         36    241ms  5.57%  6.70ms
        window sum part 2         24    128ms  2.96%  5.35ms
      solve linear systems        12    601μs  0.01%  50.0μs
      calculate flow              12    121μs  0.00%  10.1μs
    image interpolation           12    158ms  3.66%  13.2ms
    flow interpolation            12    137ms  3.16%  11.4ms
    filter inds                   12    155μs  0.00%  12.9μs
  setup                            1    177ms  4.10%   177ms
    find edge points               1    109ms  2.51%   109ms
    hist match                     1   49.7ms  1.15%  49.7ms
------------------------------------------------------------

------------------------------------------------------------
                                              Time
                                      ----------------------
           Tot / % measured:               22.6s / 100%

Section                       ncalls     time   %tot     avg
------------------------------------------------------------
pflap                              1    22.6s   100%   22.6s
  single filter pyramid level      8    22.6s   100%   2.82s
    lap                           19    9.07s  40.1%   477ms
      filtering                   19    6.39s  28.3%   336ms
      prepare A and b             19    1.05s  4.63%  55.1ms
        window sum part 1         57    664ms  2.94%  11.7ms
        window sum part 2         38    227ms  1.00%  5.96ms
      solve linear systems        19   1000ms  4.42%  52.6ms
      calculate flow              19    134ms  0.59%  7.05ms
    prefiltering                  46    6.78s  30.0%   147ms
    smoothing                     19    3.95s  17.5%   208ms
    inpainting                    19    2.34s  10.3%   123ms
      replicating borders         19   54.7ms  0.24%  2.88ms
    image interpolation           19    340ms  1.50%  17.9ms
  setup                            1   46.7ms  0.21%  46.7ms
    hist match                     1   40.9ms  0.18%  40.9ms
------------------------------------------------------------
