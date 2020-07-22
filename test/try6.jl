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
    flow_est, source_reg = polyfilter_lap(img, imgw, display=false, timer=timer)
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
method_kwargs = Dict(:timer => timer, :point_count => 35, :spacing => 35, :flow_interpolation_method => :rbf)
flow_est, source_reg, timer, results, (estim_at_inds, inds) = test_registration_alg(sparse_lap, img, imgw, flow, [25, (51,51)], method_kwargs, timer=timer)

errf
norm_err = normalize_to_zero_one(abs.(err))
sparse_estim_flow = create_sparse_flow_from_sparse(estim_at_inds, inds, size(flow_est))
addpoints(inds, labels=map(x -> string(round(x, digits=2)), norm_err))
sparse_truth_flow = create_sparse_flow_from_full(flow, inds)
showflow(sparse_truth_flow.-sparse_estim_flow)

ncc(err, vec_len.((map(ind -> flow[ind], inds) .- estim_at_inds)))
addpoints(inds, labels=map(x -> string(round(x, digits=2)), vec_len.((map(ind -> flow[ind], inds) .- estim_at_inds))))


to = timer["reg alg: sparse lap"]["sparse lap"]["prepare A and b"]

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


sections = ["reg alg: lap", "lap", "calculate flow"];
TimerOutputs.time(get_timer(old_timer, sections))/TimerOutputs.time(get_timer(new_timer, sections))

showflow(flow_est)
showflow(flow)

## PF LAP
timer=TimerOutput("pf lap");
method_kwargs =Dict(:timer => timer, :display => false, :max_repeats => 1)
flow_est, source_reg, timer, results = test_registration_alg(polyfilter_lap, img, imgw, flow, [], method_kwargs, timer=timer);

@code_warntype polyfilter_lap(img, imgw, display = false)

showflow(flow_est)
showflow(flow)
showflow(flow - flow_est, figtitle="flow - old flow")
## SPARSE PF LAP
img, imgw, flow = gen_init();
img, imgw, flow = gen_init(:lena, :uniform, flow_args=[-1 - 4im, 15]);

timer=TimerOutput("sparse pf lap");
method_kwargs = Dict(:timer => timer, :display => false, :max_repeats => 1, :point_count => 500, :spacing => 10)
flow_est, source_reg, timer, results = test_registration_alg(sparse_pflap, img, imgw, flow, [], method_kwargs, timer=timer)


showflow(flow_est)
showflow(flow)
showflow(flow-flow_est)

imgoverlay(img, source_reg)
imgshow(img - source_reg)
imgshow(img)
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
