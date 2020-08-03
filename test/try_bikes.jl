
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
