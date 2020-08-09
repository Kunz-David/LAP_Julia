using LAP_julia: normalize_to_zero_one, create_sparse_flow_from_sparse, create_sparse_flow_from_full, vec_len
using TimerOutputs

function gen_flow(img_size)

    a, b, c = randn(), randn(), randn()
    f(z) = a + b*z + c*z^2 # + d*z^3 + e*z^4

    X = ones(img_size[1]) * collect(range(0,1,length=img_size[2]))'
    Y = collect(range(0,1,length=img_size[1])) * ones(img_size[2])'
    A = Y .* im + X

    B = f.(A)
    max_len = maximum(LAP_julia.vec_len.(B))
    B = B./max_len
    return B
end

B = gen_flow((12,12)); showflow(B)

println(B)

f(1 * im)

imgshow(rand(12,123))


using TimerOutputs, LAP_julia


img, imgw, flow = gen_init(:chess, :quad, flow_args=[20])
img, imgw, flow = gen_init(:lena, :quad, flow_args=[20])

TimerOutputs.enable_debug_timings(LAP_julia)
timer = TimerOutput("Registration");
@timeit timer "polyfilter lap" begin
    flow_est, source_reg = pflap(img, imgw, display=false, timer=timer)
end
print_timer(timer)


showflow(flow_est)
showflow(flow-flow_est)

b = zeros(123,123)
@benchmark a = zeros(123,123)
@benchmark a = similar(b)
@benchmark a = Array{Float64}(undef, 123, 123)
@benchmark a = Array{Float32}(undef, 123, 123)

similar(b, (234, 234))

## fix vis methods

img = gen_chess()
img = gen_lena()
flow = ones(size(img)) .* 20 .+ ones(size(img)) .* 10 .* im
imgw = warp_img(img, -real(flow), -imag(flow))
showflow(flow)
imgshowflow(imgw, flow)
imgshow(img)

##
img, imgw, flow = gen_init()

showflow(flow)
imgshow(imgw)

## speedtest jl

img, imgw, flow = gen_init(:lena, :uniform, flow_args=[-1 - 4im, 15]);
## SPARSE LAP
timer=TimerOutput("reg alg: sparse lap");
method_kwargs = Dict(:timer => timer, :point_count => 200, :spacing => 15)
# method_kwargs = Dict(:timer => timer, :point_count => 35, :spacing => 35, :flow_interpolation_method => :rbf)
flow_est, source_reg, timer, results, (estim_at_inds, inds) = test_registration_alg(sparse_lap, img, imgw, flow, [25, (51,51)], method_kwargs, timer=timer)

norm_err = normalize_to_zero_one(abs.(err))
sparse_estim_flow = create_sparse_flow_from_sparse(estim_at_inds, inds, size(flow_est))
addpoints(inds, labels=map(x -> string(round(x, digits=2)), norm_err))
sparse_truth_flow = create_sparse_flow_from_full(flow, inds)
showflow(sparse_truth_flow.-sparse_estim_flow)

ncc(err, vec_len.((map(ind -> flow[ind], inds) .- estim_at_inds)))
addpoints(inds, labels=map(x -> string(round(x, digits=2)), vec_len.((map(ind -> flow[ind], inds) .- estim_at_inds))))


to = timer["reg alg: sparse lap"]["sparse lap"]["filtering"]

TimerOutputs.tottime(to)/10e8

showflow(flow_est-flow)

showflow(flow_est)
showflow(flow)
imgoverlay(img, source_reg)
imgoverlay(img, imgw)

## LAP
img, imgw, flow = gen_init(:lena, :tiled, flow_args=[15, 30])
timer=TimerOutput("reg alg: lap");
flow_est, source_reg, timer, results = test_registration_alg(lap, img, imgw, flow, [25, (51,51)], Dict(:timer => timer), timer=timer);


img, imgw = gen_anhir()
flow_est, source_reg = sparse_pflap(img, imgw)
showflow(flow_est)
imgshow(source_reg)
imgshow(img)
imgshow(imgw)

sections = ["reg alg: lap", "lap", "calculate flow"];
TimerOutputs.time(get_timer(old_timer, sections))/TimerOutputs.time(get_timer(new_timer, sections))

showflow(flow_est)
showflow(flow)

## PF LAP
timer=TimerOutput("pf lap");
method_kwargs =Dict(:timer => timer, :display => false, :max_repeats => 1)
flow_est, source_reg, timer, results = test_registration_alg(pflap, img, imgw, flow, method_kwargs=method_kwargs, timer=timer);

@code_warntype pflap(img, imgw, display = false)

showflow(flow_est)
showflow(flow)
showflow(flow - flow_est, figtitle="flow - old flow")
## SPARSE PF LAP
img, imgw, flow = gen_init();
img, imgw, flow = gen_init(:spaghetti, :uniform, flow_args=[-1 - 4im, 15]);

timer=TimerOutput("sparse pf lap");
method_kwargs = Dict(:timer => timer, :display => false, :max_repeats => 1, :point_count => 500, :spacing => 10)
flow_est, source_reg, timer, results = test_registration_alg(sparse_pflap, img, imgw, flow, method_kwargs=method_kwargs, timer=timer)

showflow(flow_est)
showflow(flow)
showflow(flow-flow_est)

imgoverlay(img, source_reg)
imgshow(img - source_reg)
imgshow(img)
imgshow(imgw)
imgshow(source_reg)


timer["sparse pf lap"]["single filter pyramid level"]

window_sizes = cat(collect(1:51), collect(61:10:101), collect(141:40:381), dims=1)
limit = 70
map(x -> x <= limit ? x : 0, window_sizes)
img_size = 50
img = gen_init(:chess, chess_args=[25, img_size/25])








##
filter_size = 2*filter_half_size + 1
filter_half_size = 2
k = (-filter_half_size:filter_half_size)

basis = collect(reshape(1:length(k)^2, filter_size, filter_size))

basis * k

basis' * k


##
filter_half_size = 15
filter_size = 2*filter_half_size + 1
k = (-filter_half_size:filter_half_size)

# Get the displacement vector field from the filters

basis = similar(img, (filter_size, filter_size, filter_num))
forward_ker, backward_ker = prepare_gaussian_filters(filter_half_size)
for k in 1:3
    basis[:, :, k] .= parent(broadcast(*, forward_ker[k]...))
end

basis

image_size = size(img)
all_coeffs = rand((image_size[1]^2), 3)


u1_top = zeros(image_size);
u1_bot = zeros(image_size);
u2_top = zeros(image_size);
u2_bot = zeros(image_size);

for n in 1:filter_num
    @views u1_top[:] .= u1_top[:] .- sum(transpose(basis[:, :, n]) * k) .* all_coeffs[:, n];
    @views u1_bot[:] .= u1_bot[:] .+ sum(basis[:, :, n]) .* all_coeffs[:, n];

    @views u2_top[:] .= u2_top[:] .- sum(basis[:, :, n] * k) .* all_coeffs[:, n];
    @views u2_bot[:] .= u2_bot[:] .+ sum(basis[:, :, n]) .* all_coeffs[:, n];
end

u_est = 2 .* ((im .* u1_top ./ u1_bot) .+ (u2_top ./ u2_bot));


sum(transpose(basis[:, :, 2]) * k)
sum(basis[:, :, 1] * k)
sum(basis[:, :, 2])

## use the fact that some of the sums are 0

u1_top = zeros(image_size);
u1_bot = zeros(image_size);
u2_top = zeros(image_size);
u2_bot = zeros(image_size);

for n in 1:filter_num
    if n == 1
        @views u1_top[:] .= u1_top[:] .- 0 .* all_coeffs[:, n];
        @views u1_bot[:] .= u1_bot[:] .+ 1 .* all_coeffs[:, n];

        @views u2_top[:] .= u2_top[:] .- 0 .* all_coeffs[:, n];
        @views u2_bot[:] .= u2_bot[:] .+ 1 .* all_coeffs[:, n];
    elseif n == 2
        @views u1_top[:] .= u1_top[:] .- 0 .* all_coeffs[:, n];
        @views u1_bot[:] .= u1_bot[:] .+ 0 .* all_coeffs[:, n];

        @views u2_top[:] .= u2_top[:] .- 1 .* all_coeffs[:, n];
        @views u2_bot[:] .= u2_bot[:] .+ 0 .* all_coeffs[:, n];
    elseif n == 3
        @views u1_top[:] .= u1_top[:] .- 1 .* all_coeffs[:, n];
        @views u1_bot[:] .= u1_bot[:] .+ 0 .* all_coeffs[:, n];

        @views u2_top[:] .= u2_top[:] .- 0 .* all_coeffs[:, n];
        @views u2_bot[:] .= u2_bot[:] .+ 0 .* all_coeffs[:, n];
    end
end

u_est_2 = 2 .* ((im .* u1_top ./ u1_bot) .+ (u2_top ./ u2_bot));

n = 3
sum(transpose(basis[:, :, n]) * k)
sum(basis[:, :, n])
sum(basis[:, :, n] * k)
sum(basis[:, :, n])

isapprox(u_est, u_est_2)

showflow(u_est - u_est_2)




## coeffs
B = ones(123, 3)

B[[1, 3, 12,12, 14], :] .= 12

dif = ones(123,123,3)

fhs = 2
shift = CartesianIndex(fhs, fhs)

ind = CartesianIndex(5,5)
dif[ind-shift:ind+shift, 2] .= 2

dif[ind, :]


@benchmark Dict()


function foo(img1, dict=Dict(:sone => img1))
    println(typeof(img1))
    println(typeof(dict))
end

foo(img)
img1 = img

imgfilt = zeros(size(img))
σ = 25
gaus = gaussian(σ, 2 * σ + 1)

@benchmark out = imfilter(img1, kernelfactors((gaus, gaus)), "symmetric")

@benchmark imfilter!(img1, img1, kernelfactors((gaus, gaus)), "symmetric")

out == img1
out == imgfilt


using LAP_julia: gen_uniform_flow

uni_flow = gen_uniform_flow((200, 200), 1+1im, 15)


## fair comparison matlab vs this

img, imgw, flow = gen_init(:lena, :uniform, flow_args=[-1 - 4im, 15]);

timer=TimerOutput("reg alg: lap");
flow_est, source_reg, timer, results = test_registration_alg(lap, img, imgw, flow, [25, (51,51)], Dict(:timer => timer), timer=timer);


div_coef = zeros(12,12)
div_coef[3:4,3:7] .= 1
mask = falses(12,12)
mask[2:6, 3:4] .= true
@benchmark newly_inpainted = ((div_coef .!= 0) .& (mask .== 0))
@benchmark newly_inpainted = ((div_coef .!= 0) .* (mask .== 0))


mask
div_coef

source = imgw
target = img

source_hist = adjust_histogram(source, Matching(targetimg = target))

@benchmark out = adjust_histogram($source, Matching(targetimg = $target))
@benchmark adjust_histogram!($source, Matching(targetimg = $target)) # takes 6ms on average

@benchmark (target, source) = LAP_julia.pad_images($target, $source)

@benchmark level_count = floor(Int64, log2(minimum(size(target))/8)+1)+1

@benchmark mask =

using LAP_julia: max_displacement

function foo()
    for k in 1:120
        println("iter: ", k)
        img, imgw = gen_anhir()
        try
            flow_est, source_reg = sparse_pflap(img, imgw)
        catch e
            println(e)
            return img, imgw
        end
        if max_displacement(flow_est) > 200
            return img, imgw
        end
    end
    println("all fine")
end

img, imgw = foo()

img, imgw = gen_anhir()


flow = gen_tiled_flow(size(img))
imgw = warp_img(img, real(flow), imag(flow))


imgshow(img, origin_bot=true, figtitle="target")
imgshow(imgw, origin_bot=true, figtitle="source")

max_repeats = 3;
flow_est, source_reg, figs = sparse_pflap(img, imgw, display=true, max_repeats=max_repeats)

using LAP_julia
assess_flow_quality(flow_est, flow)


@manipulate for k in slider(0:max_repeats*size(figs,1), value=1, label="k")
    if k == 0
        imgshow(imgw, origin_bot=true, figtitle="source")
    else
        reshape(permuteddimsview(figs, [2, 1, 3]), (:, 4))[k, 2]
    end
end

reshape(figs, (:, 5))

figs[1,1,2]


imgshow(source_reg, origin_bot=true, figtitle="source_reg")
imgshow(img, origin_bot=true, figtitle="target")


for k in 1:120
    try
        flow_est, source_reg, figs = sparse_pflap(img, imgw, display=true)
    catch e
        println(e)
        break
    end
end


img = gen_chess(10,10)
inds = find_edge_points(img, point_count=150, spacing=8)

imgshow(img); addpoints(inds)

##

a = zeros(Float16, 123,123);
c = ones(Float16, 5,5);
@benchmark imfilter(a, c)
b = zeros(Float64, 123,123);
d = ones(Float64, 5,5);
@benchmark imfilter(b, d)



##









img, imgw, flow = gen_init()

showflow(flow)

small_img = imresize(img, ratio = 0.5)
small_imgw = imresize(imgw, ratio = 0.5)

## small version
@btime small_flow_est, small_source_reg, rest = sparse_lap(small_img, small_imgw, 15, point_count = 500, spacing = 5)
imgshow(small_source_reg)
showflow(small_flow_est)

@btime enlarged_small_flow = imresize(small_flow_est, ratio = 2) .* 2

showflow(enlarged_small_flow.-flow)




## normal version
@btime flow_est, source_reg, rest = sparse_lap(img, imgw, 30)
imgshow(source_reg)
showflow(flow_est)

showflow(flow_est.-flow)



## classic lap
flow_est_lap, source_reg_lap = lap(img, imgw, 30)


showflow(flow_est_lap.-flow)




mutated_flow = imresize(small_flow, ratio = 3/2)


showflow(mutated_flow)
showflow(flow)
showflow(flow.-mutated_flow)



##
using CSV, FileIO

base_path = "/Users/MrTrololord/Documents/anhir/";

loc_table = CSV.read(joinpath(base_path, "location_table.csv"))
train_rows = loc_table[loc_table[:status] .== "training", :]

random_row = train_rows[rand(1:size(train_rows, 1)), :]
target_path = random_row[Symbol("Target image")]

img = load(joinpath(base_path, "dataset", target_path))
imgg = Gray.(img)

imfilter(img, ones(13,13))

eltype(img)
eltype(imgg)

sizeof(eltype(img))
sizeof(eltype(imgg))



##
target = img

level_count = floor(Int64, log2(minimum((1024,1024))/8)+1)+1

minimum(size(target))/1024

half_size_pyramid = Int64.(2 .^ range(level_count-1, stop=0, length=level_count))

1024/8

128/16



## optable

macro optable(expr)
    if expr.args[1] == :(=>) && expr.args[2] isa Int
        n = expr.args[2]
        nexpr = expr.args[3]
        name = string(nexpr.args[1])
        descr = string(nexpr)
        :(optable($(esc(nexpr)), $name, $descr, $n))
    elseif expr.args[1] == :(=>) && expr.args[2] isa String
        name = expr.args[2]
        nexpr = expr.args[3]
        descr = string(nexpr)
        :(optable($(esc(nexpr)), $name, $descr))
    else
        name = string(expr.args[1])
        descr = string(expr)
        :(optable($(esc(expr)), $name, $descr))
    end
end

function optable(op, name, descr)
    fname = joinpath("..", "assets", string(name, ".png"))
    i = 2
    while isfile(fname)
        fname = joinpath("..", "assets", string(name, i, ".png"))
        i = i + 1
    end
    out = augment(pattern, op)
    save(fname, out)
    header = length(descr) < 20 ? "Output for `$descr`" : "`$descr`"
    tbl = string(
        "Input | $header\n",
        "------|--------\n",
        "![input](../assets/testpattern.png) | ![output]($fname)\n"
    )
    Markdown.parse(tbl)
end

function optable(img1, img2, basename1, basename2, descr1, descr2)
    # save fig 1 as next in line
    imgshow(img1); i = 1;
    fname1 = joinpath("..", "..", @__DIR__, "assets", string(basename1, i, ".png"))
    while isfile(fname1)
        i = i + 1
        fname1 = joinpath("..", "..", @__DIR__, "assets", string(basename1, i, ".png"))
    end
    savefig(fname1)
    # save fig 2 as next in line
    imgshow(img2); j = 1;
    fname2 = joinpath("..", "..", @__DIR__, "assets", string(basename2, j, ".png"))
    while isfile(fname2)
        i = i + 1
        fname2 = joinpath("..", "..", @__DIR__, "assets", string(basename2, j, ".png"))
    end
    savefig(fname2)
    tbl = string(
        "$descr1 | $descr2\n",
        "------|--------\n",
        "![input]($fname1) | ![output]($fname2)\n"
    )
    Markdown.parse(tbl)
end



optable(img, imgw, "img", "imgw", "target", "source")


function homografy_flow(flow_size, max_magnitude = 10)

    a, b, c = randn(), randn(), randn()
    α, β, γ = (randn(), randn(), randn()) .* im .+ (randn(), randn(), randn())

    u(x,y)=(α*x + β*y + γ)/(a*x + b*y + c)

    X = ones(flow_size[1]) * collect(range(0,1,length=flow_size[2]))'
    Y = collect(range(0,1,length=flow_size[1])) * ones(flow_size[2])'
    A = Y .* im + X

    B = u.(X,Y)

    # max_len = maximum(LAP_julia.vec_len.(B))
    # B = B .* (max_magnitude/max_len)

end

# (Mat_<float>(3,3) << (1-rng.uniform(-0.05f, 0.05f)),
# (rng.uniform(-0.03f, 0.03f)), (rng.uniform(10.f, 20.f)),
# (rng.uniform(-0.03f, 0.03f)), (1-rng.uniform(-0.05f, 0.05f)),(rng.uniform(10.f, 20.f)),
# (rng.uniform(0.0001f, 0.0003f)), (rng.uniform(0.0001f, 0.0003f)), 1.f);


flow = homografy_flow(size(img))
showflow(flow[215:230, 100:110], skip_count = 0)


mean(flow)
median(real.(flow))
median(imag.(flow))
maximum(real.(abs.(flow)))
maximum(imag.(abs.(flow)))

findall(real.(abs.(flow)) .== maximum(real.(abs.(flow))))

imgshow(imag.(flow), origin_bot=true)




flow1 = homografy_flow(size(img))
showflow(flow1)

imgw = warp_img(img, flow)

imgshow(img)
imgshow(imgw)


timer=TimerOutput("sparse pf lap");
method_kwargs = Dict(:timer => timer, :display => false, :max_repeats => 1, :point_count => 500, :spacing => 10)
flow_est, source_reg, timer, results = test_registration_alg(sparse_pflap, img, imgw, flow1, method_kwargs=method_kwargs, timer=timer)

showflow(flow_est)

figs[1][1,1,4]

imgshow(img)



showflow(flow1)
showflow(flow_est-flow1)

timer=TimerOutput("pf lap");
method_kwargs =Dict(:timer => timer, :display => false, :max_repeats => 1)
flow_est, source_reg, timer, results = test_registration_alg(pflap, img, imgw, flow1, method_kwargs=method_kwargs, timer=timer);
