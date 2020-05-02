

point_count = 50
flow = gen_rand_flow((200, 200), 30, 40)
inds = gen_rand_points(flow, point_count, "Semigridded")


### scatter

points = reshape([ind[i] for i=1:length(inds[1]) for ind=inds], length(inds), length(inds[1]))'
samples = [flow[I] for I in inds]

meshgrid(x, y) = [repeat(x, outer=length(y)) repeat(y, inner=length(x))]
meshgrid(x::Real, y::Real) = [repeat(1:x, outer=y) repeat(1:y, inner=x)]

function test_scatter(flow_size, points, samples; method)
    gridPoints = meshgrid(flow_size...)'
    itp = interpolate(method, points, samples);
    interpolated = ScatteredInterpolation.evaluate(itp, gridPoints)
    gridded = reshape(interpolated, flow_size)
    return gridded
end

bench = @benchmark test_scatter(size(flow), points, samples; method = Polyharmonic(2))


### Dier

x = [ind[1] for ind in inds]
y = [ind[2] for ind in inds]
z1 = [real(flow[ind]) for ind in inds]
z2 = [imag(flow[ind]) for ind in inds]

function test_dier(flow_size, x, y, z1, z2)
    storage = 1000
    spl1 = Spline2D(x, y, z1; w=ones(length(x)), kx=1, ky=1, s=storage)
    real1 = evalgrid(spl1, 1:flow_size[1], 1:flow_size[2])

    spl2 = Spline2D(x, y, z2; w=ones(length(x)), kx=1, ky=1, s=storage)
    imag1 = evalgrid(spl2, 1:flow_size[1], 1:flow_size[2])

    gridded = real1 .+ imag1 .* im
    return gridded
end


@benchmark test_dier(size(flow), x, y, z1, z2)





### generate indices

"""
    gen_rand_points(flow, point_count, mode; max_skew_const=4)

Generate `point_count` random indices over `flow` depending on selection `mode`.
"""
function gen_rand_points(flow, point_count, mode; max_skew_const=4)
    if mode == "Random_points"
        ret_inds = rand(CartesianIndices(flow), point_count)
    elseif mode == "Gridded"
        points_in_dim = round(Int64, sqrt(point_count))
        point_spacing = floor.(size(flow) ./ points_in_dim)

        ret_inds = Array{CartesianIndex{2}, 1}(undef, 0)
        for ind in CartesianIndices(flow)
            if all((Tuple(ind)) .% point_spacing .== 0)
                add_ind = CartesianIndex(Tuple(ind) .- floor.(Int64, point_spacing ./ 2))
                push!(ret_inds, add_ind)
            end
        end
    elseif mode == "Semigridded"
        points_in_dim = round(Int64, sqrt(point_count))
        point_spacing = floor.(size(flow) ./ points_in_dim)

        ret_inds = Array{CartesianIndex{2}, 1}(undef, 0)
        for ind in CartesianIndices(flow)
            if all((Tuple(ind)) .% point_spacing .== 0)
                skew = (rand((-floor(Int64, point_spacing[1]/max_skew_const):floor(Int64, point_spacing[1]/max_skew_const))),
                        rand((-floor(Int64, point_spacing[1]/max_skew_const):floor(Int64, point_spacing[2]/max_skew_const))))
                add_ind = CartesianIndex(Tuple(ind) .- floor.(Int64, point_spacing ./ 2) .+ skew)
                @assert add_ind in CartesianIndices(flow)

                push!(ret_inds, add_ind)
            end
        end
    end
    return ret_inds
end



## add to df
index = 1

push!(df, Dict(:Index => index
               ,:Method => method
               ,:Mode => point_mode
               ,:Point_count => point_count
               ,:Median_speed => median(bench).time
               ,:Magnitude_of_MSE => vec_len(mse(method_flow, flow))
               ,:Angle_RMS => angle_rms(method_flow, flow)
               ,:Angle_mean => angle_mean(method_flow, flow)
               ,:Benchmark => bench
               ,:Truth_flow => showflow(flow, figtitle="Truth flow", ret="PyObject")
               ,:Method_flow => showflow(method_flow, figtitle=(method * " flow"), ret="PyObject")
               ,:Scatter_points => scatter(x, y, marker=:x)
               ))

bench = @benchmark test_scatter($flow_size, $points, $samples; method = $Polyharmonic(2))

typeof(bench)
## run


# parameters:
test_count = 3
flow_size = (200, 200)


Methods = ["Dier",
           "Scatter-Thin_plate_spline",
           "Scatter-Multiquadratic"]

Modes = ["Random_points",
         "Semigridded",
         "Gridded"]

Point_counts = [9,
                16,
                25,
                49,
                81,
                100,
                144]
import LAP_julia: helpers.mse, helpers.vec_len, helpers.angle_rms, helpers.angle_mean

let
    index = 0

    for test in test_count
        flow = gen_rand_flow(flow_size, 30, 40)
        for point_count in Point_counts
            for point_mode in Modes
                inds = gen_rand_points(flow, point_count, point_mode)
                for method in Methods
                        x = [ind[1] for ind in inds]
                        y = [ind[2] for ind in inds]
                    if method == "Dier"
                        # set inputs
                        z1 = [real(flow[ind]) for ind in inds]
                        z2 = [imag(flow[ind]) for ind in inds]
                        # do a benchmark
                        bench = @benchmark test_dier($x, $y, $z1, $z2)
                        # eval with current input
                        method_flow = test_dier(x, y, z1, z2)
                    elseif method == "Scatter-Thin_plate_spline"
                        # set inputs
                        points = reshape([ind[i] for i=1:length(inds[1]) for ind=inds], length(inds), length(inds[1]))'
                        samples = [flow[I] for I in inds]
                        # do a benchmark
                        bench = @benchmark test_scatter($flow_size, $points, $samples; method = $Polyharmonic(2))
                        # eval with current input
                        method_flow = test_scatter(flow_size, points, samples; method = Polyharmonic(2))
                    elseif method == "Scatter-Multiquadratic"
                        # set inputs
                        points = reshape([ind[i] for i=1:length(inds[1]) for ind=inds], length(inds), length(inds[1]))'
                        samples = [flow[I] for I in inds]
                        # do a benchmark
                        bench = @benchmark test_scatter($flow_size, $points, $samples; method = $Multiquadratic(2))
                        # eval with current input
                        method_flow = test_scatter(flow_size, points, samples; method = Multiquadratic(2))
                    end
                    # ADD TO DATA FRAME
                    index = index + 1
                    push!(df, Dict(:Index => index,
                                   :Method => method,
                                   :Mode => point_mode,
                                   :Point_count => point_count,
                                   :Median_speed => median(bench).time,
                                   :Magnitude_of_MSE => vec_len(mse(method_flow, flow)),
                                   :Angle_RMS => angle_rms(method_flow, flow),
                                   :Angle_mean => angle_mean(method_flow, flow),
                                   :Benchmark => bench,
                                   :Truth_flow => showflow(flow, figtitle="Truth flow", ret="PyObject"),
                                   :Method_flow => showflow(method_flow, figtitle=(method * " flow"), ret="PyObject"),
                                   :Scatter_points => scatter(x, y, marker=:x)))
                    println("at index", index)
                end
            end
        end
    end
end # let
