
using LinearAlgebra

function e_med(u, u_est)
    return median(abs.(LAP_julia.vec_len.(u .- u_est))[:])
end

function e_mean(u, u_est)
    return mean(abs.(LAP_julia.vec_len.(u .- u_est))[:])
end

## init

# lena
using TestImages
img = testimage("lena_gray")
img = Float32.(img)
wanted_max_displacement = 15
constant_flow = false

imgshow(img)
save("../PolyFilter_LAP/target.png", img)


##
# create 3 flows
if constant_flow == false
    # one
    base1 = [1 + im 2 - im;  2 - im 1 + im]
    tile_size1 = round(Int64, size(img, 1)/2)
    base_flow1 = repeat(base1, inner=(tile_size1, tile_size1))
    flow1 = LAP_julia.smooth_with_gaussian(base_flow1, 40)
    mag_ratio1 = wanted_max_displacement/LAP_julia.max_displacement(flow1)
    flow1 = flow1 .* mag_ratio1
    showflow(flow1)

    # two
    base2 = [1 + im 2 - 2im;  2 - im 1]
    tile_size2 = round(Int64, size(img, 1)/2)
    base_flow2 = repeat(base2, inner=(tile_size2, tile_size2))
    flow2 = LAP_julia.smooth_with_gaussian(base_flow2, 50)
    mag_ratio2 = wanted_max_displacement/LAP_julia.max_displacement(flow2)
    flow2 = flow2 .* mag_ratio2
    showflow(flow2)

    # three
    base3 = [-1 + im 2 + 2im;  2 - im im]
    tile_size3 = round(Int64, size(img, 1)/2)
    base_flow3 = repeat(base3, inner=(tile_size3, tile_size3))
    flow3 = LAP_julia.smooth_with_gaussian(base_flow3, 70)
    mag_ratio3 = wanted_max_displacement/LAP_julia.max_displacement(flow3)
    flow3 = flow3 .* mag_ratio3
    showflow(flow3)
else
     # one
    flow1 = repeat([-1 + 1im], inner=(size(img)))
    mag_ratio1 = wanted_max_displacement/LAP_julia.max_displacement(flow1);
    flow1 = flow1 .* mag_ratio1;

    # two
    flow2 = repeat([-12 + 1im], inner=(size(img)))
    mag_ratio2 = wanted_max_displacement/LAP_julia.max_displacement(flow2);
    flow2 = flow2 .* mag_ratio2;

    # one
    flow3 = repeat([-1 - 4im], inner=(size(img)))
    mag_ratio3 = wanted_max_displacement/LAP_julia.max_displacement(flow3);
    flow3 = flow3 .* mag_ratio3;
end


## test lap methods


#results
pflap_times = zeros(3)
pflap_e_meds = zeros(3)
pflap_e_means = zeros(3)

# one
flow = flow1; num = 1;
imgw = warp_img(img, -real(flow), -imag(flow));
u_est, source_reg = pflap(img, imgw, display=false);
bench = @benchmark pflap(img, imgw, display=false)
pflap_times[num] = median(bench.times) / 10^9
pflap_e_meds[num] = e_med(flow, u_est)
pflap_e_means[num] = e_mean(flow, u_est)

# two
flow = flow2; num = 2;
imgw = warp_img(img, -real(flow), -imag(flow));
u_est, source_reg = pflap(img, imgw, display=false);
bench = @benchmark pflap(img, imgw, display=false)
pflap_times[num] = median(bench.times) / 10^9
pflap_e_meds[num] = e_med(flow, u_est)
pflap_e_means[num] = e_mean(flow, u_est)


# three
flow = flow3; num = 3;
imgw = warp_img(img, -real(flow), -imag(flow));
u_est, source_reg = pflap(img, imgw, display=false);
bench = @benchmark pflap(img, imgw, display=false)
pflap_times[num] = median(bench.times) / 10^9
pflap_e_meds[num] = e_med(flow, u_est)
pflap_e_means[num] = e_mean(flow, u_est)
showflow(u_est)

pflap_times
pflap_e_meds
pflap_e_means

## single_lap

# specifics
fhs = 23
window = [47, 47]

#results
lap_times = zeros(3)
lap_e_meds = zeros(3)
lap_e_means = zeros(3)

# one
flow = flow1; num = 1;
imgw = warp_img(img, -real(flow), -imag(flow));
u_est = lap(img, imgw, fhs, window)
bench = @benchmark lap(img, imgw, fhs, window)
lap_times[num] = median(bench.times) / 10^9
lap_e_meds[num] = e_med(flow, u_est)
lap_e_means[num] = e_mean(flow, u_est)

# two
flow = flow2; num = 2;
imgw = warp_img(img, -real(flow), -imag(flow));
u_est = lap(img, imgw, fhs, window)
bench = @benchmark lap(img, imgw, fhs, window)
lap_times[num] = median(bench.times) / 10^9
lap_e_meds[num] = e_med(flow, u_est)
lap_e_means[num] = e_mean(flow, u_est)

# three
flow = flow3; num = 3;
imgw = warp_img(img, -real(flow), -imag(flow));
u_est = lap(img, imgw, fhs, window)
bench = @benchmark lap(img, imgw, fhs, window)
lap_times[num] = median(bench.times) / 10^9
lap_e_meds[num] = e_med(flow, u_est)
lap_e_means[num] = e_mean(flow, u_est)

lap_times
lap_e_meds
lap_e_means

## sparse single lap

# specifics
fhs = 23
window = [47, 47]

#results
splap_times = zeros(3)
splap_e_meds = zeros(3)
splap_e_means = zeros(3)

# one
flow = flow1; num = 1;
imgw = warp_img(img, -real(flow), -imag(flow));
u_est = new_alg(img, imgw, fhs, window)
bench = @benchmark new_alg(img, imgw, fhs, window)
splap_times[num] = median(bench.times) / 10^9
splap_e_meds[num] = e_med(flow, u_est)
splap_e_means[num] = e_mean(flow, u_est)

# two
flow = flow2; num = 2;
imgw = warp_img(img, -real(flow), -imag(flow));
u_est = new_alg(img, imgw, fhs, window)
bench = @benchmark new_alg(img, imgw, fhs, window)
splap_times[num] = median(bench.times) / 10^9
splap_e_meds[num] = e_med(flow, u_est)
splap_e_means[num] = e_mean(flow, u_est)

# three
flow = flow3; num = 3;
imgw = warp_img(img, -real(flow), -imag(flow));
u_est = new_alg(img, imgw, fhs, window)
bench = @benchmark new_alg(img, imgw, fhs, window)
splap_times[num] = median(bench.times) / 10^9
splap_e_meds[num] = e_med(flow, u_est)
splap_e_means[num] = e_mean(flow, u_est)

splap_times
splap_e_meds
splap_e_means

# classic
function lap(img, imgw, fhs, window_size)
    classic_estim = single_lap(img, imgw, fhs, window_size)
    LAP_julia.inpaint_nans!(classic_estim)
    LAP_julia.smooth_with_gaussian(classic_estim, window_size)
    return classic_estim
end

# new
function new_alg(img, imgw, fhs, window_size)
    mask = parent(padarray(trues(size(img).-(2*fhs, 2*fhs)), Fill(false, (fhs, fhs), (fhs, fhs))))
    inds = find_edge_points(img, mask=mask)
    points = LAP_julia.inds_to_points(inds)
    new_estim = single_lap_at_points(img, imgw, fhs, window_size, points, 3)
    full_new_estim = interpolate_flow(new_estim, inds)
    return full_new_estim
end


## sparse pf lap

#results
sppflap_times = zeros(3)
sppflap_e_meds = zeros(3)
sppflap_e_means = zeros(3)

# one
flow = flow1; num = 1;
imgw = warp_img(img, -real(flow), -imag(flow));
u_est, source_reg = sparse_pflap(img, imgw, display=false)
bench = @benchmark new_alg(img, imgw, fhs, window)
sppflap_times[num] = median(bench.times) / 10^9
sppflap_e_meds[num] = e_med(flow, u_est)
sppflap_e_means[num] = e_mean(flow, u_est)

# two
flow = flow2; num = 2;
imgw = warp_img(img, -real(flow), -imag(flow));
u_est, source_reg = sparse_pflap(img, imgw, display=false)
bench = @benchmark new_alg(img, imgw, fhs, window)
sppflap_times[num] = median(bench.times) / 10^9
sppflap_e_meds[num] = e_med(flow, u_est)
sppflap_e_means[num] = e_mean(flow, u_est)

# three
flow = flow3; num = 3;
imgw = warp_img(img, -real(flow), -imag(flow));
u_est, source_reg = sparse_pflap(img, imgw, display=false)
bench = @benchmark new_alg(img, imgw, fhs, window)
sppflap_times[num] = median(bench.times) / 10^9
sppflap_e_meds[num] = e_med(flow, u_est)
sppflap_e_means[num] = e_mean(flow, u_est)
showflow(u_est)
showflow(flow .- u_est)

imgshow(img .- source_reg)
imgoverlay(img, source_reg)



sppflap_times
sppflap_e_meds
sppflap_e_means

## make comparation table:

lap = [lap_times, lap_e_meds, lap_e_means]
pflap = [pflap_times, pflap_e_meds, pflap_e_means]
splap = [splap_times, splap_e_meds, splap_e_means]
sppflap = [sppflap_times, sppflap_e_meds, sppflap_e_means]

columns = ["Time", "e_med", "e_mean"]
lap_avg = round.(map(mean, lap), digits=3)
pflap_avg = round.(map(mean, pflap), digits=3)
splap_avg = round.(map(mean, splap), digits=3)
sppflap_avg = round.(map(mean, sppflap), digits=3)

println(columns)
map(println, [pflap_avg, lap_avg, sppflap_avg, splap_avg])



## results
map(println, [lap_avg, pflap_avg, splap_avg, sppflap_avg])
"LENA"
["Time", "e_med", "e_mean"]
[1.198, 1.384, 1.806]
[4.668, 0.852, 1.012]
[0.232, 1.331, 1.835]
[0.205, 0.546, 0.619]

"const"
map(println, [pflap_avg, lap_avg, sppflap_avg, splap_avg])
["Time", "e_med", "e_mean"]
6.664 & 0.116 & 0.721
1.019 & 2.553 & 3.149
0.194 & 0.428 & 0.988
0.196 & 2.285 & 2.524
