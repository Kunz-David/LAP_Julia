
using LAP_julia, TimerOutputs, DataFrames, TableView, JLD2, FileIO
##
sp_lap_point_count = DataFrame(
    index = Int[],
    diag_size = Int[],
    whs = Int[],
    timer = TimerOutput[],
    median_results = Dict{String,Float64}[],
    point_count = Int[],
    spacing = Int[],
    points_found = Float16[]
    );

window_half_sizes = cat(collect(15:15:51), dims=1)
diag_sizes = [400, 600, 900, 1200];
point_counts = vcat(collect(50:150:700), collect(800:300:1550))
spacings = vcat(5, 15, collect(35:15:40), 50, collect(60:20:155))

##
using LAP_julia: fun_on_dict_values, resize_to_diag_size
using BenchmarkTools

img_count = 5;
errors = Array{Any}(undef, 0)
imgs = Array{Image}(undef, 5)
# prepare images:
for k in 1:img_count
    imgs[k], _ = gen_anhir(mutate=false)
end

imgshow(imgs[2])
imgshow(resize_to_diag_size(imgs[2], 450))
size(resize_to_diag_size(imgs[2], 500))

# erase old_results
df = sp_lap_point_count;

let index = 0
    for diag_size in diag_sizes
        for whs in window_half_sizes
            for point_count in point_counts, spacing in spacings
                timer = TimerOutput("reg alg: sp lap")
                results_dicts = Array{Dict}(undef, img_count)
                inds_array = Array{Array{Any}}(undef, img_count)

                for k in 1:img_count
                    img = resize_to_diag_size(imgs[k], diag_size)
                    whs_limited = round(Int, min(whs, min(size(img)...)/4))
                    rand_flow = gen_quad_flow(size(img), whs_limited)
                    imgw = warp_img(img, real(rand_flow), imag(rand_flow))
                    try
                        _, _, timer, results_dicts[k], (_, inds_array[k]) = test_registration_alg(sparse_lap, img, imgw, rand_flow, [whs_limited], Dict(:timer => timer, :point_count => point_count, :spacing => spacing), timer=timer, display=false)
                        # @assert 1==0
                    catch e
                        println(e)
                        push!(errors, [e, img, imgw, rand_flow, whs_limited])
                        continue
                    end
                end
                # if no data gathered:
                if any(x-> 0 .== x, length.([results_dicts, inds_array]))
                    println("SKIPPED: diag_size:", diag_size, " whs: ", whs)
                    continue
                end

                median_results = fun_on_dict_values(filter_defined(results_dicts), median)
                avg_points_found = LAP_julia.mean(map(x -> length(x), filter_defined(inds_array)))

                index = index + 1
                println("at index: ", index, " diag_size:", diag_size, " whs: ", whs)
                push!(df, Dict(:index => index,
                               :diag_size => diag_size,
                               :whs => whs,
                               :timer => timer,
                               :median_results => median_results,
                               :point_count => point_count,
                               :spacing => spacing,
                               :points_found => avg_points_found
                               ))
            end # point_count in point_counts, spacing in spacings
        end
    end
end

defined_perms(arraylike) = filter(k -> isassigned(arraylike, k), 1:length(arraylike))
filter_defined(arraylike) = arraylike[defined_perms(arraylike)]

results_dicts[defined_perms(results_dicts)]
filter_defined(results_dicts)

errors

results_dicts = Array{Dict}(undef, img_count)

results_dicts[2] = Dict("sd" => 12)

filter(k -> isassigned(results_dicts, k), 1:length(results_dicts))

isassigned(results_dicts, 2)

filter(, results_dicts)

# df

# win_sum1_df = df[(df.reg_fun .== :sparse_lap_win_sum1), :];
# win_sum3_df = df[(df.reg_fun .== :sparse_lap), :]; # [1, :timer]
#
# TimerOutputs.ncalls(df[(df.reg_fun .== :sparse_lap), :][1, :timer][sections[1]][sections[2]][sections[3]])
#
# get_timer(win_sum3_df[1, :timer], sections)

import LAP_julia: sparse_lap


win_sum3_df[1, :timer]


sections = ["reg alg: sp lap", "sparse lap", "prepare A and b"]
sections = ["reg alg: sp lap", "sparse lap", "filtering"]

fig, axs = subplots(6, sharex=true, figsize=(10,10))
for (diag_size, ax) in zip(diag_sizes[1:end-1], axs)
    ax.set_title("Image Size: " * string(diag_size))

    n = count(k -> k < (diag_size/4), window_half_sizes)
    timers3 = df[(df.reg_fun .== :sparse_lap) .& (df.diag_size .== diag_size), :][:, :timer][1:n]

    y3 = map(x -> TimerOutputs.time(get_timer(x, sections))/
                  (TimerOutputs.ncalls(get_timer(x, sections)) * 10e8), timers3)
    x = filter(x -> x < (diag_size/4), window_half_sizes)[1:n]

    println(length(y3), "  ", length(x))

    xscale("log")
    yscale("log")
    subplots_adjust(top=1.5)
    ax.plot(x, y3)
    xlabel("whs [pixels]")
    ax.set_ylabel("time [s]")
end
gcf()



showtable(df)

TimerOutputs.time(get_timer(win_sum3_df[:, :timer][1], ["reg alg: sp lap", "sparse lap", "prepare A and b"]))
TimerOutputs.tottime(get_timer(win_sum3_df[:, :timer][1], ["reg alg: sp lap", "sparse lap", "prepare A and b"]))
TimerOutputs.tottime(get_timer(win_sum3_df[:, :timer][1], ["reg alg: sp lap", "interpolate flow"]))
TimerOutputs.time(win_sum3_df[:, :timer][1]["reg alg: sp lap"]["interpolate flow"])


print_timer(win_sum3_df[:, :timer][1])


function get_timer(timer, sections)
    private_sections = copy(sections)
    try
        timer[private_sections[1]]
    catch e
        println(e)
        return TimerOutput()
    end
    out = timer
    while private_sections != []
        out = out[popfirst!(private_sections)]
    end
    return out
end

full_flow_estim, source_reg, flow_estim_at_inds, inds = sparse_lap_filt_one_by_one(img, imgw, 12)

inds

using ImageFiltering: KernelFactors.gaussian
using LAP_julia: inds_to_points, filt_onebyone!, window_sum3_at_inds!, multi_mat_div_using_qr
function sparse_lap_filt_one_by_one(target,
                                    source,
                                    fhs,
                                    window_size=(2fhs+1, 2fhs+1);
                                    spacing::Integer=15,
                                    point_count::Integer=100,
                                    timer::TimerOutput=TimerOutput("sparse_lap"),
                                    display::Bool=false,
                                    flow_interpolation_method::Symbol=:quad,
                                    base_method_kwargs=Dict(:timer => timer, :display => display))

    mask = parent(padarray(trues(size(target).-(2*fhs, 2*fhs)), Fill(false, (fhs, fhs), (fhs, fhs))))
    @timeit_debug timer "find edge points" begin
        inds = find_edge_points(target, spacing=spacing, number=point_count, mask=mask)
    end
    @timeit_debug timer "sparse lap" begin
        flow_estim_at_inds, inds = single_lap_at_points_filt_one_by_one(target, source, fhs, window_size, inds; base_method_kwargs...)
    end
    if all(isnan, flow_estim_at_inds)
        @timeit_debug timer "interpolate flow" begin
            full_flow_estim = zeros(size(target)) .* im .+ zeros(size(target))
        end
    else
        @timeit_debug timer "interpolate flow" begin
            full_flow_estim = interpolate_flow(flow_estim_at_inds, inds, size(target), method=flow_interpolation_method)
        end
    end
    @timeit_debug timer "generate source_reg" begin
        source_reg = warp_img(source, real(full_flow_estim), imag(full_flow_estim))
    end
    return full_flow_estim, source_reg, flow_estim_at_inds, inds
end

function single_lap_at_points_filt_one_by_one(target::Image,
                                              source::Image,
                                              filter_half_size::Integer,
                                              window,
                                              inds::Array{CartesianIndex,1};
                                              filter_count::Integer=3,
                                              timer::TimerOutput=TimerOutput("Sparse LAP"),
                                              display::Bool=false)

    image_size = size(target)
    pixel_count = length(target)
    ind_count = length(inds)
    filter_size = 2 * filter_half_size + 1

    # Prepare filters:
    # Calculate separable filters from basis:
    sigma = (filter_half_size + 2) / 4
    centered_inds = centered(-filter_half_size:filter_half_size)
    gaus = gaussian(sigma, filter_size)
    gaus_der = gaus .* centered_inds .* (-1)/sigma^2
    gaus_der_flip = reverse(gaus_der, dims=1)

    # prepare 2D filter basis for later use
    basis = similar(target, (filter_size, filter_size, filter_count))

    # dif = (points of forward image - points of backward image) for each filter
    dif = similar(target, (image_size..., filter_count))
    points = inds_to_points(inds)

    @timeit_debug timer "filtering" begin
        if filter_count == 3
            # temporary place to store filtered images
            # tmp_filtered_1 = fill(NaN, image_size)
            # tmp_filtered_2 = fill(NaN, image_size)
            tmp_filtered_1 = similar(target)
            tmp_filtered_2 = similar(target)

            # basis 1 - (gaus, gaus) both forward and backward
            kernf_1 = kernelfactors((gaus, gaus))
            basis[:, :, 1] = broadcast(*, kernf_1...)
            filt_onebyone!(tmp_filtered_1, target, kernf_1, filter_half_size, points)
            filt_onebyone!(tmp_filtered_2, source, kernf_1, filter_half_size, points)
            # imfilter!(tmp_filtered_1, target, kernf_1, "symmetric")
            # imfilter!(tmp_filtered_2, source, kernf_1, "symmetric")
            dif[:, :, 1] = tmp_filtered_2 - tmp_filtered_1

            # basis 2 - (gaus, gaus_der_flip) as forward, (gaus, gaus_der) as backward
            kernf_2f = kernelfactors((gaus, gaus_der_flip))
            basis[:, :, 2] = broadcast(*, kernf_2f...)
            filt_onebyone!(tmp_filtered_1, target, kernf_2f, filter_half_size, points)
            # imfilter!(tmp_filtered_1, target, kernf_2f, "symmetric")
            kernf_2b = kernelfactors((gaus, gaus_der))
            filt_onebyone!(tmp_filtered_2, source, kernf_2b, filter_half_size, points)
            # imfilter!(tmp_filtered_2, source, kernf_2b, "symmetric")
            dif[:, :, 2] = tmp_filtered_2 - tmp_filtered_1

            # basis 3 - (gaus_der_flip, gaus) as forward, (gaus_der, gaus) as backward
            kernf_3f = kernelfactors((gaus_der_flip, gaus))
            basis[:, :, 3] = broadcast(*, kernf_3f...)
            filt_onebyone!(tmp_filtered_1, target, kernf_3f, filter_half_size, points)
            # imfilter!(tmp_filtered_1, target, kernf_3f, "symmetric")
            kernf_3b = kernelfactors((gaus_der, gaus))
            filt_onebyone!(tmp_filtered_2, source, kernf_3b, filter_half_size, points)
            # imfilter!(tmp_filtered_2, source, kernf_3b, "symmetric")
            dif[:, :, 3] = tmp_filtered_2 - tmp_filtered_1
        end
    end # "filtering"

    dif = reshape(dif, (:, filter_count))

    # prepare matrices for linear system of equations
    A = similar(target, (filter_count-1, filter_count-1, ind_count))
    b = similar(target, (filter_count-1, ind_count))

    @timeit_debug timer "prepare A and b" begin
        # TODO: check whether window_sum_around_points calculates only around points and not on whole window.
        for k in 1:filter_count-1
            for l in k:filter_count-1
                @timeit_debug timer "window sum part 1" begin
                    # @views window_sum_around_points!(A[k, l, :], dif[:, k+1] .* dif[:, l+1], image_size, window, points)
                    @views window_sum3_at_inds!(A[k, l, :], dif[:, k+1] .* dif[:, l+1], image_size, window, inds)
                    A[l, k, :] = A[k, l, :]
                end
            end
            @timeit_debug timer "window sum part 2" begin
                # @views window_sum_around_points!(b[k, :], dif[:, k+1] .* dif[:, 1] .* (-1), image_size, window, points)
                @views window_sum3_at_inds!(b[k, :], dif[:, k+1] .* dif[:, 1] .* (-1), image_size, window, inds)
            end
        end
    end # "prepare A and b"

    # Perform Gauss elimination on all pixels in parallel:
    # coeffs will be of shape: pixel_count, filter_count-1
    @timeit_debug timer "multi mat div" begin
        @views coeffs = multi_mat_div_using_qr(A, b)
    end

    # adding ones so that all base filters have their coefficients even the first one
    all_coeffs = [ones(ind_count) coeffs]

    k = (-filter_half_size:filter_half_size)

    # Get the displacement vector field from the filters
    @timeit_debug timer "calculate flow" begin
        u_est_at_inds = similar(target, ComplexF64, ind_count)
        @views u_est_at_inds[:] = 2 .* ((im .* (.-1 .* all_coeffs[:, 3]) ./ all_coeffs[:, 1]) .+ ((.-1 .* all_coeffs[:, 2]) ./ all_coeffs[:, 1]));
    end # "calculate flow"

    # dont use estimations whose displacement is larger than the filter_half_size
    displacement_mask = (real(u_est_at_inds).^2 .+ imag(u_est_at_inds).^2) .<= filter_half_size^2

    return u_est_at_inds[displacement_mask], inds[displacement_mask]

end
