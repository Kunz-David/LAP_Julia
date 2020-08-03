
using LAP_julia

viewimg(img) = colorview(Gray, img)


img_size = (250,255)



view(img)
view(imgw)

timer = TimerOutput(string(args_dict["reg_alg"]))
method_kwargs = Dict(:display => true, :timer => timer)
img, imgw, flow = gen_init(:spaghetti, img_args=[img_size], flow_args=[25])
flow_est, source_reg, timer, time_in_secs = test_registration_alg(sparse_pflap,
                                                                  img,
                                                                  imgw,
                                                                  flow,
                                                                  timer=method_kwargs[:timer],
                                                                  method_kwargs=method_kwargs)

showflow(flow_est)
showflow(flow.*(-1))

count(isnan.([real.(flow_est) imag.(flow_est)]))


LAP_julia.vec_len.(flow_est)
view(source_reg)
view(img)
view(imgw)



##
base_path_bikes = "/Users/MrTrololord/Downloads/bikes"

img = load(joinpath(base_path_bikes, "img1.png"))
for k in 2:2
    imgw = load(joinpath(base_path_bikes, "img$(k).png"))
    global imgw
end

imgw

img = Float64.(Gray.(img))
imgw = Float64.(Gray.(imgw))

timer = TimerOutput("sparse_pflap");
method_kwargs = Dict(:display => false, :timer => timer)
flow_est, source_reg, timer, time_in_secs = time_reg_alg(pflap, img, imgw, timer=method_kwargs[:timer], method_kwargs=method_kwargs)

view(imgw)
view(source_reg)
view(img)

showflow(flow_est)

##










base_path_zhang = "/Users/MrTrololord/Downloads/dataset/"
fnames = ["yard", "PGH", "testchart"]
imgs = Array{Image}(undef, size(fnames))
imgws = Array{Image}(undef, size(fnames))


for (k, name) in zip(1:length(fnames), fnames)
    imgs[k] = load(joinpath(base_path_zhang, "$(name)_t.jpg"))
    imgws[k] = load(joinpath(base_path_zhang, "$(name)_s.jpg"))
    global img, imgw
end

img = imgs[3]
imgw = imgws[3]

img = Float64.(Gray.(img))
imgw = Float64.(Gray.(imgw))

##
timer = TimerOutput("pflap");
method_kwargs = Dict(:display => false, :timer => timer)
flow_est, source_reg, timer, time_in_secs = time_reg_alg(pflap, img, imgw, timer=method_kwargs[:timer], method_kwargs=method_kwargs)

source_reg_zeros = warp_img(imgw, -real.(flow_est), -imag.(flow_est), border_strat=:zeros)
view(source_reg_zeros)

view(imgw)
view(source_reg)
viewimg(img)

showflow(imresize(flow_est, (200,200)))


##
# @info "size from $(size)"
(img, ratio), (imgw, ratio) = resize_to_diag_size(img, 400), resize_to_diag_size(imgw, 400)

img, imgw = Gray.(img), Gray.(imgw)

timer = TimerOutput("sparse_pflap");
method_kwargs = Dict(:display => false, :timer => timer)
flow_est, source_reg, timer, time_in_secs = time_reg_alg(sparse_pflap_psnr, imgw, img, timer=method_kwargs[:timer], method_kwargs=method_kwargs)

source_reg_zeros = warp_img(img, -real.(flow_est), -imag.(flow_est), border_strat=:zeros)
viewimg(source_reg_zeros)
viewimg(imgw)

flow_est

##
# load float from text

function parse_numbers(s)
    pieces = split(s, ' ', keepempty=false)
    map(pieces) do piece
        parse(Float64, piece)
    end
end

function load_H(path)
    lines = readlines(path)
    H = Array{Float64}(undef, (3,3))
    for (k, line) in zip(1:3, lines)
        H[k,:] = parse_numbers(line)
    end
    return H
end

function make_flow_from_H(H, size)
    warp_perspective(y, x) = (H[1,1]*x + H[1,2]*y + H[1,3])/(H[3,1]*x + H[3,2]*y + H[3,3]),
                             (H[2,1]*x + H[2,2]*y + H[2,3])/(H[3,1]*x + H[3,2]*y + H[3,3])
    flow = Array{Complex{Float64},2}(undef, size)
    X = ones(size[1]) * collect(range(0,1,length=size[2]))'
    Y = collect(range(0,1,length=size[1])) * ones(size[2])'
    for ind in CartesianIndices(flow)
        calc = [Y[ind], X[ind]] .- warp_perspective(Y[ind], X[ind])
        flow[ind] = calc[1] + im*calc[2]
    end
    return flow
end

H = load_H(joinpath(base_path_bikes, "H1to2p"))
out = make_flow_from_H(H, size(img))

showflow(out)

img
imgw

asdf = warp_img(img, -imag.(out), -real.(out), border_strat=:zeros)

Gray.(img) .- Gray.(asdf)

##
using CoordinateTransformations, StaticArrays, ImageTransformations, LinearAlgebra;
struct Homography{T} <: AbstractAffineMap
    m::SMatrix{3, 3, T, 9}
end

# having T in the RHS forces calling the constructor again
Homography(m::AbstractMatrix{}) = Homography(m);
function (trans::Homography{M})(x::SVector{3}) where M <: Real
    out = trans.m * x;
    out = out / out[end];
    SVector{2}(out[1:2])
end

function (trans::Homography{M})(x::SVector{2}) where M <: Real
    trans(SVector{3}([x[1], x[2], 1.0]))
end

function (trans::Homography{M})(x::CartesianIndex{2}) where M <: Real
    trans(SVector{3}([collect(x.I); 1]))
end

function (trans::Homography{M})(x::Tuple{Int, Int}) where M <: Real
    trans(CartesianIndex{2}(x))
end

function (trans::Homography{M})(x::Array{CartesianIndex{2}, 1}) where M <: Real
    CartesianIndex{2}.([tuple(y...) for y in trunc.(Int, collect.(trans.(x)))])
end

function Base.inv(trans::Homography)
    i = inv(trans.m);
    Homography(i ./ i[end])
end

Homo = Homography{Float64}(H)

img1_warp = ImageTransformations.warp(img, Homo)

imgw



##

## posun homografii

base_path_bikes = "/Users/MrTrololord/Downloads/bikes"

img = load(joinpath(base_path_bikes, "img1.png"))
k = 2
imgw = load(joinpath(base_path_bikes, "img$(k).png"))

function parse_numbers(s)
    pieces = split(s, ' ', keepempty=false)
    map(pieces) do piece
        parse(Float64, piece)
    end
end

function load_H(path)
    lines = readlines(path)
    H = Array{Float64}(undef, (3,3))
    for (k, line) in zip(1:3, lines)
        H[k,:] = parse_numbers(line)
    end
    return H
end

function make_flow_from_H(H, size)
    warp_perspective(y, x) = (H[1,1]*x + H[1,2]*y + H[1,3])/(H[3,1]*x + H[3,2]*y + H[3,3]),
                             (H[2,1]*x + H[2,2]*y + H[2,3])/(H[3,1]*x + H[3,2]*y + H[3,3])
    flow = Array{Complex{Float64},2}(undef, size)
    X = ones(size[1]) * collect(range(0,1,length=size[2]))'
    Y = collect(range(0,1,length=size[1])) * ones(size[2])'
    for ind in CartesianIndices(flow)
        calc = [Y[ind], X[ind]] .- warp_perspective(Y[ind], X[ind])
        flow[ind] = calc[1] + im*calc[2]
    end
    return flow
end

H = load_H(joinpath(base_path_bikes, "H1to2p"))
out = make_flow_from_H(H, size(img))

warped_img = warp_img(img, -imag.(out), -real.(out), border_strat=:zeros)
imgw

## posun mym alg

viewimg(img) = colorview(Gray, img)

img, imgw = Gray.(img), Gray.(imgw)

timer = TimerOutput("sparse_pflap");
method_kwargs = Dict(:display => false, :timer => timer)
flow_est, source_reg, timer, time_in_secs = time_reg_alg(sparse_pflap_psnr, imgw, img, timer=method_kwargs[:timer], method_kwargs=method_kwargs)

source_reg_zeros = warp_img(img, -real.(flow_est), -imag.(flow_est), border_strat=:zeros)
viewimg(source_reg_zeros)
viewimg(imgw)


# posun mym alg je ma lepsi vysledky?? nesmysl -> neco spatne
