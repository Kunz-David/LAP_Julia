
## posun homografii

base_path_bikes = "/Users/MrTrololord/Downloads/ox_affine/images/leuven/"

img = load(joinpath(base_path_bikes, "img1.png"))
k = 3
imgw = load(joinpath(base_path_bikes, "img$(k).png"))

H = load_H(joinpath(base_path_bikes, "H1to$(k)p"))
out = make_flow_from_H(H, size(img))


img, imgw_large = LAP_julia.pad_images(img, imgw)
warped_img = warp_img(img, real.(out), imag.(out), border_strat=:zeros)
# viewimg(warped_img)
imgw_large

# show warp:
showflow(imresize(out, ratio=0.4))

## posun mym alg

viewimg(img) = colorview(Gray, img)

img, imgw = Gray.(img), Gray.(imgw)
# img, imgw = imresize(img, (287, 410)), imresize(imgw, (287, 410))


method = sparse_pflap
timer = TimerOutput(string(method));
method_kwargs = Dict(:display => false, :timer => timer, :match_source_histogram => true)
flow_est, source_reg, timer, time_in_secs = time_reg_alg(method, imgw, img, timer=method_kwargs[:timer], method_kwargs=method_kwargs)

source_reg_zeros = warp_img(img, -real.(flow_est), -imag.(flow_est), border_strat=:zeros)
viewimg(source_reg_zeros)
imgw_large

# ground  truth
warped_img

showflow(imresize(flow_est.*(-1), ratio=0.4))
# showflow(imresize(view(flow_est.*(-1), 1:50, 1:50), ratio=0.4))

viewimg(source_reg_zeros.-imgw)


imgoverlay(source_reg_zeros, imgw, figtitle="prekryv s moji registraci")
imgoverlay(warped_img, imgw, figtitle="prekryv s posunem z homografie")

# posun mym alg je ma lepsi vysledky?? nesmysl -> neco spatne, neumim posunout homografii.. :(
