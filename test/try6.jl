
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

timer=TimerOutput("reg alg: sparse lap");
flow_est, source_reg, timer, results = test_registration_alg(sparse_lap, img, imgw, flow, [12, (25,25)], Dict(:timer => timer), timer=timer)

to = timer["reg alg: sparse lap"]["sparse lap"]["prepare A and b"]

TimerOutputs.tottime(to)/10e8

showflow(flow_est)
showflow(flow)

timer=TimerOutput("reg alg: lap");
test_registration_alg(lap, img, imgw, flow, [25, (51,51)], Dict(:timer => timer), timer=timer)


window_sizes = cat(collect(1:51), collect(61:10:101), collect(141:40:381), dims=1)
limit = 70
map(x -> x <= limit ? x : 0, window_sizes)
img_size = 50
img = gen_init(:chess, chess_args=[25, img_size/25])
