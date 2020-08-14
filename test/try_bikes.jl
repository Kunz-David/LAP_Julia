
## posun homografii

base_path_bikes = "/Users/MrTrololord/Downloads/ox_affine/images/bikes/"
# base_path_bikes = "/BIRL/mount/ox_affine/images/leuven/"

img = load(joinpath(base_path_bikes, "img1.png"))
k = 4
imgw = load(joinpath(base_path_bikes, "img$(k).png"))

H = LAP_julia.load_H(joinpath(base_path_bikes, "H1to$(k)p"))
out = LAP_julia.make_flow_from_H(H, size(img))

flip_out = reverse(out, dims = 1)
flip_out = real.(flip_out) .- imag(flip_out) .* 1im

showflow(imresize(flip_out, ratio=0.1))

img, imgw_large = LAP_julia.pad_images(img, imgw)
warped_img = warp_img(img, real.(out), imag.(out), border_strat=:zeros)
# viewimg(warped_img)
imgw_large

# show warp:
maximum(LAP_julia.vec_len.(out))
showflow(imresize(out, ratio=0.1))
showflow(imresize(reverse(out, dims = 1).*(-1), ratio=0.1))
showflow(imresize(reverse(out, dims = 1), ratio=0.1))

img
imgw

#reverse
imgww = warp_img(imgw, -real.(out), -imag.(out), border_strat=:zeros)
img

## posun mym alg

viewimg(img) = colorview(Gray, img)

img, imgw = Gray.(img), Gray.(imgw)
img, imgw = Float64.(Gray.(img)), Float64.(Gray.(imgw))
# img, imgw = imresize(img, (287, 410)), imresize(imgw, (287, 410))


method = sparse_pflap_psnr
timer = TimerOutput(string(method));
method_kwargs = Dict(:display => false, :timer => timer, :match_source_histogram => true, :point_count => 200)
flow_est, source_reg, timer, time_in_secs = test_registration_alg(method, img, imgw, out.*(-1), timer=method_kwargs[:timer], method_kwargs=method_kwargs);

method = pflap
timer = TimerOutput(string(method));
method_kwargs = Dict(:display => false, :timer => timer, :match_source_histogram => true, :prefilter => true)
flow_est, source_reg, timer, time_in_secs = test_registration_alg(method, img, imgw, out.*(-1), timer=method_kwargs[:timer], method_kwargs=method_kwargs);



source_reg_zeros = warp_img(imgw, -real.(flow_est), -imag.(flow_est), border_strat=:zeros)
img
viewimg(source_reg_zeros)
imgw_large

# with landmarks
imgshow(source_reg_zeros, origin_bot=true)
addpoints(transformed_source_inds)

imgshow(imgw, origin_bot=true)
addpoints(source_inds)

# ground  truth
warped_img

showflow(imresize(flow_est, ratio=0.1))
showflow(imresize(out, ratio=0.1))
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




### get max displacement

ox_folder = "/Users/MrTrololord/Downloads/ox_affine/images/"

folders = ["bikes", "leuven", "trees", "ubc"]
# folder = "trees"
for folder in folders
    k = 2
    H = LAP_julia.load_H(joinpath(ox_folder, folder, "H1to$(k)p"))
    out = LAP_julia.make_flow_from_H(H, size(img))
    cur_max = maximum(LAP_julia.vec_len.(out))
    @info cur_max
    max_disp = cur_max
    min_max_disp = cur_max

    for k in 3:6
        H = LAP_julia.load_H(joinpath(ox_folder, folder, "H1to$(k)p"))
        out = LAP_julia.make_flow_from_H(H, size(img))
        cur_max = maximum(LAP_julia.vec_len.(out))
        @info cur_max
        max_disp = max(cur_max, max_disp)
        min_max_disp = min(cur_max, min_max_disp)
    end

    @info folder max_disp min_max_disp
end

### get example images

ox_folder = "/Users/MrTrololord/Downloads/ox_affine/images/"
plots_folder = "/Users/MrTrololord/Google_Drive/cvut/bakalarka/plots"

folders = ["bikes", "leuven", "trees", "ubc"]
method = sparse_pflap

for folder in folders
    k = 3
    img = load(joinpath(ox_folder, folder, "img1.png"))
    imgw = load(joinpath(ox_folder, folder, "img$(k).png"))
    H = LAP_julia.load_H(joinpath(ox_folder, folder, "H1to$(k)p"))
    ground_flow = LAP_julia.make_flow_from_H(H, size(img))
    img_save_path = joinpath(plots_folder, "$(folder)_target.png")
    save(img_save_path, img)
    imgw_save_path = joinpath(plots_folder, "$(folder)_source.png")
    save(imgw_save_path, imgw)

    flip_out = reverse(ground_flow, dims = 1)
    flip_out = real.(flip_out) .- imag(flip_out) .* 1im

    showflow(flip_out, figtitle="Ground Truth Deformation")
    ground_flow_save_path = joinpath(plots_folder, "$(folder)_ground_truth_flow.png")
    savefig(ground_flow_save_path)

    timer = TimerOutput(string(method));
    method_kwargs = Dict(:display => false, :timer => timer, :match_source_histogram => true)
    img, imgw = Gray.(img), Gray.(imgw)
    flow_est, source_reg, timer, time_in_secs = test_registration_alg(method, img, imgw, ground_flow.*(-1), timer=method_kwargs[:timer], method_kwargs=method_kwargs);

    flip_flow_est = reverse(flow_est, dims = 1)
    flip_flow_est = real.(flip_flow_est) .- imag(flip_flow_est) .* 1im

    showflow(flip_flow_est, figtitle="Estimated Deformation")
    ground_flow_save_path = joinpath(plots_folder, "$(folder)_estimated_flow_with_$(string(method)).png")
    savefig(ground_flow_save_path)

    @info folder
end






------------------------------------------------------------
                                              Time
                                      ----------------------
           Tot / % measured:               18.1s / 97.8%

Section                       ncalls     time   %tot     avg
------------------------------------------------------------
sparse_pflap_psnr                  1    17.7s   100%   17.7s
  single filter pyramid level      8    17.5s  98.7%   2.18s
    sparse lap                    15    16.5s  93.4%   1.10s
      filtering                   15    15.8s  89.0%   1.05s
      prepare A and b             15    768ms  4.34%  51.2ms
        window sum part 1         45    473ms  2.67%  10.5ms
        window sum part 2         30    295ms  1.66%  9.82ms
      solve linear systems        15    483μs  0.00%  32.2μs
      calculate flow              15   94.3μs  0.00%  6.29μs
    image interpolation           15    612ms  3.46%  40.8ms
    flow interpolation            15    274ms  1.55%  18.3ms
    filter inds                   15    283μs  0.00%  18.8μs
  setup                            1    236ms  1.33%   236ms
    find edge points               1    153ms  0.86%   153ms
    hist match                     1   65.1ms  0.37%  65.1ms
------------------------------------------------------------


------------------------------------------------------------
                                              Time
                                      ----------------------
           Tot / % measured:               44.1s / 100%

Section                       ncalls     time   %tot     avg
------------------------------------------------------------
pflap                              1    43.9s   100%   43.9s
  single filter pyramid level      8    43.8s   100%   5.47s
    lap                           16    19.2s  43.7%   1.20s
      filtering                   16    13.9s  31.6%   866ms
      prepare A and b             16    2.00s  4.56%   125ms
        window sum part 1         48    1.07s  2.44%  22.3ms
        window sum part 2         32    491ms  1.12%  15.4ms
      solve linear systems        16    1.75s  3.98%   109ms
      calculate flow              16    330ms  0.75%  20.6ms
    prefiltering                  40    14.0s  31.9%   350ms
    smoothing                     16    7.98s  18.2%   499ms
    inpainting                    16    1.92s  4.37%   120ms
      replicating borders         16    293ms  0.67%  18.3ms
    image interpolation           16    493ms  1.12%  30.8ms
  setup                            1   90.9ms  0.21%  90.9ms
    hist match                     1   82.7ms  0.19%  82.7ms
------------------------------------------------------------
