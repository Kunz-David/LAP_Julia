## initialize

using BenchmarkTools



###

img1, imgw1, flow1 = gen_init(flow_args=(20,60));
flow_est1, source_reg1, figs1 = sparse_pflap(img1, imgw1)


showflow(flow_est1, figtitle="estim")
showflow(flow1, figtitle="orig")

### porovnani o kolik je single_lap_at_points rychlejsi nez single_lap podle velikosti posunu


img, imgw, flow = gen_init(flow_args=(4,100));
showflow(flow, figtitle="orig")

# vars
fhs = 5
window_size = [11,11]

# classic
clasic_estim = single_lap(img, imgw, fhs, window_size)
showflow(clasic_estim, figtitle="clasic_estim")

# new
spacing = 35
point_count = 30
mask = parent(padarray(trues(size(target).-(2*fhs, 2*fhs)), Fill(false, (fhs, fhs), (fhs, fhs))))
inds = find_edge_points(target, spacing=spacing, number=point_count, mask=mask)
points = LAP_julia.inds_to_points(inds)
new_estim = single_lap_at_points(img, imgw, fhs, window_size, 3, points)
showflow(new_estim, disp_type=:sparse, figtitle="new_estim")

full_new_estim = interpolate_flow(new_estim, inds)
showflow(full_new_estim, figtitle="full_new_estim")
addpoints(inds)

imgshow(img); addpoints(inds)



# classic
function lap(img, imgw, fhs, window_size)
    classic_estim = single_lap(img, imgw, fhs, window_size)
    LAP_julia.inpaint_nans!(classic_estim)
    LAP_julia.smooth_with_gaussian(classic_estim, window_size)
    return classic_estim
end

# new
function new_alg(img, imgw, fhs, window_size)
    spacing = 35
    point_count = 35
    mask = parent(padarray(trues(size(img).-(2*fhs, 2*fhs)), Fill(false, (fhs, fhs), (fhs, fhs))))
    inds = find_edge_points(img, spacing=spacing, number=point_count, mask=mask)
    points = LAP_julia.inds_to_points(inds)
    new_estim = single_lap_at_points(img, imgw, fhs, window_size, 3, points)
    full_new_estim = interpolate_flow(new_estim, inds)
    return full_new_estim
end


function get_erroring_params_from_new(maxd)
    # vars
    runs = 0
    while true
        img, imgw, flow = gen_init(flow_args=(maxd, 100));
        fhs = maxd+2
        window_size = [fhs*2+1,fhs*2+1]

        output = new_alg(img, imgw, fhs, window_size)
        if length(output) < 200
            return output
        end
        runs = 1 + runs
        println(runs)
    end
end

out = get_erroring_params_from_new(5)

function compare_speed_dep_on_max_displacement(maxd)

    # vars
    img, imgw, flow = gen_init(flow_args=(maxd, 100));
    fhs = maxd+2
    window_size = [fhs*2+1,fhs*2+1]

    new_flow = new_alg(img, imgw, fhs, window_size)
    classic_flow = lap(img, imgw, fhs, window_size)

    new_bench = @benchmark $new_alg($img, $imgw, $fhs, $window_size)
    classic_bench = @benchmark $lap($img, $imgw, $fhs, $window_size)

    classic_speed = median(classic_bench.times)
    classic_mse = mse(classic_flow, flow)

    new_speed = median(new_bench.times)
    new_mse = mse(new_flow, flow)

    new_mse_pix = mse(warp_img(imgw, -real(new_flow), -imag(new_flow)), flow)
    classic_mse_pix = mse(warp_img(imgw, -real(classic_flow), -imag(classic_flow)), flow)

    println("CLASSIC: speed: ", classic_speed, " mse: ", classic_mse, "mse pix :", classic_mse_pix)
    println("NEW: speed: ", new_speed, " mse: ", new_mse, "mse pix :", new_mse_pix)


    println("mse pix quility increase: ", classic_mse_pix/new_mse_pix)

    speedup = classic_speed / new_speed
    quality_compare = classic_mse / new_mse

    println("speedup depending on points: ", speedup)
    println("quality increase depending on points: ", quality_compare)
    return speedup, quality_compare
end

function run_multiple(fun, run_count; args=[])
    outs = Any[]
    for k in 1:run_count
        println()
        out = fun(args...)
        push!(outs, out)
    end
    # make sure the out is a tuple or array
    if typeof(outs[1]) <: Union{Tuple, Array}
        output_size = length(outs[1])
    else
        output_size = 1
        outs = map(x -> (x), outs)
    end
    # average them
    ret = (LAP_julia.mean([outs[l][k] for l in 1:output_size]) for k in 1:output_size)
    return outs
end

using LinearAlgebra

function end_point_error(u, u_est)
    as, df = LAP_julia.mse.(u, u_est)
end
u = flow


function e_med(u, u_est)
    retrun median((LAP_julia.vec_len.(u) .- LAP_julia.vec_len.(u_est))[:])
end



(real(flow), imag(flow)) .- (real(u_est), imag(u_est)) .^2

as, df = (real(flow), imag(flow)) .- (real(u_est), imag(u_est))
as .^2

function angular_error(u, u_est)
    return 1/cos(
        (1 + dot(u[:], u_est[:])) /
        (sqrt(1 + dot(u[:], u[:])) * sqrt(1 + dot(u_est[:], u_est[:]))))
end

LAP_julia.vec_len(angular_error(flow, u_est))
end_point_error(flow, u_est)

out1, out2 = compare_speed_dep_on_max_displacement(4)



showflow(out1)

run_multiple(compare_speed_dep_on_max_displacement, 4, args=[15])


ar = [[1212, 11], [-1212, 11]]
outs = ar
output_size = 2
[LAP_julia.mean([outs[l][k] for l in 1:output_size]) for k in 1:output_size]

displacements = [2:3:38]
run







### single lap at point showcase:


# vars
maxd = 5
img, imgw, flow = gen_init(flow_args=(maxd, 100));
fhs = maxd+2
window_size = [fhs*2+1,fhs*2+1]
showflow(flow, figtitle="orig")

single_points_full_flow = new_alg(img, imgw, fhs, window_size)
showflow(single_points_full_flow, figtitle="estim using single points")

single_full_flow = lap(img, imgw, fhs, window_size)
showflow(single_full_flow, figtitle="estim using classic single")


# show dif
showflow(single_full_flow .- orig, figtitle="classic single - orig")
showflow(single_points_full_flow .- orig, figtitle="single points - orig")
