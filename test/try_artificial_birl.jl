
viewimg(img) = colorview(Gray, img)

base_path = "/Users/MrTrololord/Downloads/head_mri/"
target_images_path = joinpath(base_path, "images", "target")
k = 1

img = load(joinpath(target_images_path, "img$(k).png"))
flow = gen_quad_flow(size(img), 40)
imgw = warp_img(img, -real.(flow), -imag.(flow))



method = sparse_pflap
timer = TimerOutput(string(method));
method_kwargs = Dict(:display => false, :timer => timer, :match_source_histogram => true)
flow_est, source_reg, timer, time_in_secs = time_reg_alg(method, img, imgw, timer=method_kwargs[:timer], method_kwargs=method_kwargs);


source_reg_zeros = warp_img(imgw, -real.(flow_est), -imag.(flow_est), border_strat=:zeros)
img


## see points

imgshow(source_reg, figtitle = "Source Registered")
imgshow(img, figtitle = "Target")

showflow(flow.*(-1), figtitle = "Truth Flow")
showflow(flow_est, figtitle = "Estimated Flow")


inds, mag = find_edge_points(img, debug=true)

imgshow(mag)
addpoints(inds)
