
## posun homografii

base_path_bikes = "/Users/MrTrololord/Downloads/ox_affine/images/bikes/"
base_path_bikes = "/BIRL/mount/ox_affine/images/leuven/"

img = load(joinpath(base_path_bikes, "img1.png"))
k = 2
imgw = load(joinpath(base_path_bikes, "img$(k).png"))

H = LAP_julia.load_H(joinpath(base_path_bikes, "H1to$(k)p"))
out = LAP_julia.make_flow_from_H(H, size(img))

img, imgw_large = LAP_julia.pad_images(img, imgw)
warped_img = warp_img(img, real.(out), imag.(out), border_strat=:zeros)
# viewimg(warped_img)
imgw_large

# show warp:
showflow(imresize(out, ratio=0.4))

#reverse
imgww = warp_img(imgw, -real.(out), -imag.(out), border_strat=:zeros)
img

## posun mym alg

viewimg(img) = colorview(Gray, img)

img, imgw = Gray.(img), Gray.(imgw)
img, imgw = Float64.(Gray.(img)), Float64.(Gray.(imgw))
# img, imgw = imresize(img, (287, 410)), imresize(imgw, (287, 410))


method = sparse_pflap
method = pflap
timer = TimerOutput(string(method));
method_kwargs = Dict(:display => false, :timer => timer, :match_source_histogram => true, :prefilter => true)
method_kwargs = Dict(:display => false, :timer => timer, :match_source_histogram => true)
method_kwargs = Dict(:display => false, :timer => timer, :match_source_histogram => true, :spacing => 20)
flow_est, source_reg, timer, time_in_secs = test_registration_alg(method, img, imgw, out, timer=method_kwargs[:timer], method_kwargs=method_kwargs);

source_reg_zeros = warp_img(imgw, -real.(flow_est), -imag.(flow_est), border_strat=:zeros)
img
viewimg(source_reg_zeros)
imgw_large

# ground  truth
warped_img

showflow(imresize(flow_est.*(-1), ratio=0.4))
# showflow(imresize(view(flow_est.*(-1), 1:50, 1:50), ratio=0.4))

imgoverlay(source_reg_zeros, imgw, figtitle="prekryv s moji registraci")
imgoverlay(warped_img, imgw, figtitle="prekryv s posunem z homografie")

## see points

imgshow(source_reg, figtitle = "Source Registered")
imgshow(img, figtitle = "Target")

showflow(flow.*(-1), figtitle = "Truth Flow")
showflow(flow_est, figtitle = "Estimated Flow")


inds, mag = find_edge_points(img, debug=true, spacing = 20)
imgshow(mag)
addpoints(inds)




### test with landmarks
imgshow(img, origin_bot=true)
imgshow(imgw, origin_bot=true)
showflow(out)

transfored_imgw = warp_img(imgw, real.(out), imag.(out), border_strat=:zeros)


const COUNT = 500

target_inds = LAP_julia.gen_rand_points(img, COUNT, "Gridded")
for k in 2:6
    # get flow
    transform_path = joinpath(base_path_bikes, "H1to$(k)p")
    cur_H = LAP_julia.load_H(transform_path)
    cur_flow = LAP_julia.make_flow_from_H(cur_H, size(img))
    # shift and filter landmarks
    global target_inds = LAP_julia.get_valid_landmarks(cur_flow, target_inds)
end

target_inds

imgshow(img, origin_bot=true)
addpoints(target_inds)
addpoints(transformed_source_inds)

source_inds = LAP_julia.move_landmarks(collect(transpose(LAP_julia.inds_to_points(target_inds))), flow_est)

transformed_source_inds = LAP_julia.move_landmarks(collect(source_inds), (-1).*flow_est)

imgshow(imgw, origin_bot=true)
addpoints(source_inds)

transfored_imgw = warp_img(imgw, -real.(flow_est), -imag.(flow_est), border_strat=:zeros)
imgshow(transfored_imgw, origin_bot=true)
addpoints(transformed_source_inds)

#
# iimg, iimgw, iflow = gen_init(:lena, :uniform, flow_args=[1+1im, 20])
# imgshow(iimg, origin_bot=true)
# imgshow(iimgw, origin_bot=true)
# showflow(iflow)
#
# asdf = warp_img(iimg, -real.(iflow), -imag.(iflow), border_strat=:zeros)
# imgshow(asdf, origin_bot=true)
#
# asdf = warp_img(iimgw, real.(iflow), imag.(iflow), border_strat=:zeros)
# imgshow(asdf, origin_bot=true)
#
# imgoverlay_v2(iimg, asdf)
# imgoverlay_v2(iimg, asdf)

img
imgw
