
tile_size = 50
board_size = 4 # must be even
mini_board = [zeros(tile_size, tile_size) ones(tile_size, tile_size);
              ones(tile_size, tile_size) zeros(tile_size, tile_size)]

chessboard = repeat(mini_board, outer=(convert(Integer,(board_size/2)), convert(Integer,(board_size/2))))
img = chessboard

###
using Random, PyPlot
using ScatteredInterpolation
using Dierckx

# params:
point_count = 100


flow = gen_rand_flow(size(img), 30, 40)


rand_inds = rand(CartesianIndices(flow), point_count)

flow[rand_inds[1]]

## ScatteredInterpolation:
using BenchmarkTools, ScatteredInterpolation

points = reshape([ind[i] for i=1:length(rand_inds[1]) for ind=rand_inds], length(rand_inds), length(rand_inds[1]))'
samples = [flow[I] for I in rand_inds]

idx = rand_inds
@benchmark reduce((x,y)->hcat(x, [y.I...]), idx[2:end]; init=[idx[1].I...])



@benchmark reshape([ind[i] for i=1:length(rand_inds[1]) for ind=rand_inds], length(rand_inds), length(rand_inds[1]))'


# fun
flow_size = size(flow)
x1 = 1:flow_size[1]
y1 = 1:flow_size[2]

X = repeat(x1, flow_size[1])[:]
Y = repeat(y1', flow_size[1])[:]
gridPoints = [X Y]'

itp = interpolate(method, points, samples);
interpolated = ScatteredInterpolation.evaluate(itp, gridPoints)
gridded = reshape(interpolated, flow_size)

# fun end

# methods: RBFs = [Gaussian(2), Multiquadratic(2), InverseQuadratic(2), InverseMultiquadratic(2)]
function test_scatter(points, samples, smoothness; method)
    flow_size = size(flow)
    x1 = 1:flow_size[1]
    y1 = 1:flow_size[2]

    X = repeat(x1, flow_size[1])[:]
    Y = repeat(y1', flow_size[1])[:]
    gridPoints = [X Y]'

    itp = interpolate(method, points, samples);
    interpolated = ScatteredInterpolation.evaluate(itp, gridPoints)
    gridded = reshape(interpolated, flow_size)
    return gridded
end

@btime gridded_scatter = test_scatter($points, $samples, $0)
@benchmark gridded_scatter = test_scatter($points, $samples, $0, method=Shepard())
@benchmark gridded_scatter = test_scatter($points, $samples, $0)
gridded_scatter = test_scatter(points, samples, 0, method=Multiquadratic(2))
gridded_scatter = test_scatter(points, samples, 0, method=Shepard())

showflow(gridded_scatter, figtitle="scatter")

showflow(flow, figtitle="orig")

import LAP_julia: helpers.mse, helpers.angle_rms, helpers.angle_mean, helpers.vec_len

mse(gridded_scatter, flow)
vec_len(mse(gridded_scatter, flow))
angle_rms(gridded_scatter, flow)
angle_mean(gridded_scatter, flow)


sum(flow - gridded_scatter)/length(flow)

showflow(flow)
x = [ind[1] for ind in rand_inds]
y = [ind[2] for ind in rand_inds]
scatter(x, y, marker=:x)
gcf()


## Dierckx
using Dierckx

x = float([ind[1] for ind in rand_inds])
y = float([ind[2] for ind in rand_inds])
z1 = [real(flow[ind]) for ind in rand_inds]
z2 = [imag(flow[ind]) for ind in rand_inds]

function test_dier(x, y, z1, z2)
    storage = 1000
    spl1 = Spline2D(x, y, z1; w=ones(length(x)), kx=1, ky=1, s=storage)
    real1 = evalgrid(spl1, 1:size(flow)[1], 1:size(flow)[2])

    spl2 = Spline2D(x, y, z2; w=ones(length(x)), kx=1, ky=1, s=storage)
    imag1 = evalgrid(spl2, 1:size(flow)[1], 1:size(flow)[2])

    gridded = real1 .+ imag1 .* im
    return gridded
end

@benchmark gridded_dier = test_dier($x, $y, $z1, $z2)
gridded_dier = test_dier(x, y, z1, z2)


showflow(gridded_dier, figtitle="dier")
showflow(flow)


### IrregularInterpolate

# using PyPlot
#
# f(x, y) = (x+y-1)^2
#
# x = float([ind[1] for ind in rand_inds]) ./ 200
# y = float([ind[2] for ind in rand_inds]) ./ 200
# a = f.(x, y)
#
# p = Main.IrregularInterpolate.PolyharmonicSpline(2, [x y], a)  # 3 is the order
#
# inds = CartesianIndices(flow)
# x2 = float(vec([I[1] for I in inds])) ./ 200
# y2 = float(vec([I[2] for I in inds])) ./ 200
# data = Main.IrregularInterpolate.interpolate(p, [x2 y2])
#
# gridded = reshape(data, size(flow))
#
# imgshow(gridded)
#
#
# figure()
# plot_surface([I[1] for I in inds], [I[2] for I in inds], gridded)
# gcf()

###

using PyPlot

x = float([ind[1] for ind in rand_inds])
y = float([ind[2] for ind in rand_inds])
z1 = [real(flow[ind]) for ind in rand_inds]
z2 = [imag(flow[ind]) for ind in rand_inds]

function test_code(x, y, z1, z2)
    p1 = Main.IrregularInterpolate.PolyharmonicSpline(2, [x y], z1)  # 3 is the order
    p2 = Main.IrregularInterpolate.PolyharmonicSpline(2, [x y], z2)  # 3 is the order

    inds = CartesianIndices(flow)
    x2 = float(vec([I[1] for I in inds]))
    y2 = float(vec([I[2] for I in inds]))
    data_real = Main.IrregularInterpolate.interpolate(p1, [x2 y2])
    data_imag = Main.IrregularInterpolate.interpolate(p2, [x2 y2])

    gridded = reshape((data_real .+ data_imag .* im), size(flow))
    return gridded
end


@benchmark gridded_code = test_code($x, $y, $z1, $z2)
gridded_code = test_code(x, y, z1, z2)

showflow(gridded_code)

x = [ind[1] for ind in rand_inds]
y = [ind[2] for ind in rand_inds]
scatter(x, y, marker=:x)
gcf()

showflow(flow)

figure()
plot_surface([I[1] for I in inds], [I[2] for I in inds], imag(gridded), alpha=0.5)
gcf()


z1 = [real(point) for point in z]
z2 = [imag(point) for point in z]


### compare results

# prepare data:
flow = gen_rand_flow(size(img), 30, 60)

# set params:
point_count = 10
rand_inds = rand(CartesianIndices(flow), point_count) # random points


# prepare inputs
points = reshape([ind[i] for i=1:length(rand_inds[1]) for ind=rand_inds], length(rand_inds), length(rand_inds[1]))'
samples = [flow[I] for I in rand_inds]
x = float([ind[1] for ind in rand_inds])
y = float([ind[2] for ind in rand_inds])
z1 = [real(flow[ind]) for ind in rand_inds]
z2 = [imag(flow[ind]) for ind in rand_inds]


# scatter
@benchmark gridded_scatter = test_scatter($points, $samples, $0)
gridded_scatter = test_scatter(points, samples, 0)
# Dierckx
@benchmark gridded_dier = test_dier($x, $y, $z1, $z2)
gridded_dier = test_dier(x, y, z1, z2)
# code
b = @benchmark gridded_code = test_code($x, $y, $z1, $z2)
gridded_code = test_code(x, y, z1, z2)

# show flows:
showflow(flow, figtitle="orig")
showflow(gridded_scatter, figtitle="scatter")
showflow(gridded_dier, figtitle="dier")
showflow(gridded_code, figtitle="code")


# add points
s = scatter(x, y, marker=:x)
gcf()



###
using DataFrames

Methods = ["Dier",
           "Scatter - Thin plate spline",
           "Scatter - Multiquadratic"]

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



###
# gridded
test_flow = CartesianIndices(ones(10,10))

flow_size = size(test_flow)
point_count = 10
points_in_dim = round(Int64, sqrt(point_count))
point_spacing = floor.(flow_size ./ points_in_dim)


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
        print(size(ret_inds))
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

points = gen_rand_points(ones(100, 100), 40, "Semigridded")

imgshow(ones(100, 100))
x = [ind[1] for ind in points]
y = [ind[2] for ind in points]
scatter(x, y, marker=:x)
gcf()

# make Semigridded

###

scatter_interp_bench = DataFrame(Methods = )

df = DataFrame(A = 1:2, a = ["ss", "sdfd"], B = [bench, bench], B2 = [showflow(flow, ret="PyObject"), showflow(flow, ret="PyObject")], C =[ 1, 4])

display(df)


showflow(flow, figtitle="orig", ret="PyObject")

import PyCall.PyObject

scatter_interp_bench = DataFrame(
                                Index = Int[],
                                Method = String[],
                                Mode = String[],
                                Point_count = Int[],
                                Median_speed = Float64[],
                                Magnitude_of_MSE = Float64[],
                                Angle_RMS = Float64[],
                                Angle_mean = Float64[],
                                Benchmark = BenchmarkTools.Trial[],
                                Truth_flow = PyObject[],
                                Method_flow = PyObject[],
                                Scatter_points = PyObject[]
                                )

df = scatter_interp_bench
df = DataFrame(A = Int[], B = String[])

push!(df, (1, "shit"))


showtable(df)

insert!(df, 3, undef, :aaArea)

df.aaArea = [s, s]

# test
A = rand(50000, 30)
d = DataFrame(A)

tmp = d[d.x1 .> 0.7, :]
tmp === d[d.x1 .> 0.7, :]

arr = [12.3545235, 12.323, 2.123412]

mat = [12; 12]

round.(arr, digits=2)

vcat([3, 4, 5], arr)

## test try catch
arr = []

for k in 1:12
    try
        tmp = 12
        if k == 4
            @assert true == false
        end
    catch e
        push!(arr, (k, e))
        continue;
    end
    println(tmp)
    println(k)
end

try
    @assert true == false
catch e
    @warn e
    typeof(e)
end


function f(x::Type, y) where Type <: Real
    x .+ y
end

function g2(x::Type, y::mat) where {Type <:Real, mat <: Matrix{<:Complex{<:Real}}}
    x .+ y
end

@code_warntype f(12, fl)

g2(12, 12)

fl = gen_rand_flow((12,12), 12)

typeof(fl)


point_count_df = DataFrame(Point_c = Int[],
                           Method = String[],
                           Median_speed_mean = Float64[],
                           Magnitude_of_MSE_mean = Float64[],
                           Angle_RMS_mean = Float64[],
                           Angle_Mean_mean = Float64[])

for points in Point_counts
    for method in Methods
        tmp_df = df[((df.Point_count == points) .& (df.Method .== method)), :]
        tmp_df = select(tmp_df, :Median_speed, :Magnitude_of_MSE, :Angle_RMS, :Angle_mean)
        test_count = length(tmp_df.Median_speed)
        tmp_df = describe(tmp_df)
        retrieved_data = select(tmp_df, :mean, copycols=false)
        push!(poin_count_df, [points method retrieved_data])
    end
end

bench = @benchmark 12+1

df = DataFrame(a = Matrix[])

push!(df, Dict(:a => flow))

df[:, :a][1]

typeof(flow)

df = DataFrame()

flow = gen_rand_flow((100,100), 12)
plot = showflow(flow)

@save "save_df.JLD2" df flow
