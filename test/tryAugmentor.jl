using Augmentor

img, imgw, flow = gen_init(:chess, flow_args=[15]);
showflow(flow)

pl = FlipX() * FlipY() |> Zoom(0.9:0.1:1.2)

Augmentor.augment_impl(pl)

pl = Rotate(-10:101)
pl = ElasticDistortion(10, 10, 0.3, 4, 3, true)

img_processed = augment(img, pl); imgshow(img_processed)

img_processed = augment((img, positions), pl); imgshow(img_processed)
imgshow(img)

flow_processed = augment(flow, pl); showflow(flow_processed)

A = [1 2 3; 4 5 6; 7 8 9]

Juno.@enter augment(A, pl)

pl = Rotate(90)

B = augment(A, pl)
imgshow(A)
imgshow(B)

pl = ElasticDistortion(3, 3, 0.25, 2, 3, true, true)


img_processed = augment(img, plborder)
imgshow(img_processed)

imgshow(img)


mutable struct Pos
    x::Float64
    y::Float64
end

Base.:+(a::Pos, b::Pos) = Pos(a.x + b.x, a.y + b.y)
Base.:*(a::Any, b::Pos) = Pos(a * b.x, a * b.y)
Base.:*(b::Pos, a::Any) = Pos(a * b.x, a * b.y)
Base.zero(Pos) = Pos(0, 0)
Pos(a::Pos) = a

Base.show(io::IO, m::Pos) = print(io, "(", m.x, ", ", m.y, ")")

zero(Pos)

1/2 * Pos(1, 2) + 1/2 * Pos(1, 1)

function ind_to_pos(ind)
    return Pos(ind[1], ind[2])
end

positions = map(ind_to_pos, CartesianIndices(img))

function make_flow(positions::Array{Pos,2})
    fun = (ind, pos) -> pos.x - ind[1] + im * (pos.y - ind[2])
    foo = (ind, pos) -> ind[1] - pos.x + im * (ind[2] - pos.y)
    return map(foo, CartesianIndices(positions), positions) .* (-1)
end

inds = CartesianIndices(positions)
inds[1,1][1]
((ind, pos) -> ind[1] - pos.x + im * (ind[2] - pos.y))(CartesianIndex(1,1), Pos(2,2))

augflow = make_flow(end_positions)
showflow(augflow)
imgshow(img_processed)


pl = Rotate(9)
begin
    img_processed, end_positions = augment((img, positions), pl);
    augflow = make_flow(end_positions);
    imgshowflow(img_processed, augflow)
end

augflow
imgshowflow(img_processed, augflow)

img, imgw, flow = gen_init(:chess, flow_args=[20, 130]);
imgshowflow(imgw, flow)
addpoints([(50, 0)])

imgshow(imgw)

showflow(ones(3,3) + im * ones(3,3))

addpoints([(50, 0)])

imgw[54,1]

end_positions


## profiling
const to = TimerOutput();

reset_timer!();

timeit_debug_enabled()

@timeit "section" sleep(0.02)
@timeit "section2" sleep(0.1)
@timeit_debug "debug time" sleep(1)

TimerOutputs.enable_debug_timings(LAP_julia)

print_timer()

img, imgw, flow = gen_init()

u_est, source_reg = pflap(img, imgw, display=false);

showflow(flow)
showflow(u_est)
imgshow(source_reg)
imgshow(img)


@profiler pflap(img, imgw, display=false);

## single lap

function lap(img, imgw, fhs, window_size)
    classic_estim = single_lap(img, imgw, fhs, window_size)
    LAP_julia.inpaint_nans!(classic_estim)
    LAP_julia.smooth_with_gaussian(classic_estim, window_size)
    LAP_julia.smooth_with_gaussian(classic_estim, window_size)
    return classic_estim
end

fhs = 23
window = [47, 47]
u_est = lap(img, imgw, fhs, window);

showflow(u_est)

@time lap(img, imgw, fhs, window)

# profile
Profile.init()
Profile.init(10^10, 0.00001)
Juno.@profiler lap(img, imgw, fhs, window);


## test A
C = [1 2 3; 4 5 6; 7 8 19]
C[3, :]

k, l = size(C)

A = rand(2, 3, 4)

A[:, :, 1][1,2]

A[2,1,:]




C = [1 2 3; 4 5 6; 7 8 19]

function gem2d(A, b)
    for k in axes(A,1)
        for l in k+1:size(A, 2)
            ratio = A[l, k] / A[k, k]
            println(ratio)
            A[l, :] .-= ratio .* A[k, :]
            b[l] -= ratio * b[k]
        end
    end
    return A, b
end

gem2d(C, [1, 2, 3])



ratios = zeros(size(A, 3))

function gem3d!(A, b)
    ratios = zeros(size(A, 3))
    for k in axes(A,1)
        for l in k+1:size(A, 2)
            ratios .= A[l, k, :] ./ A[k, k, :]
            A[l, :, :] .-= ratios' .* A[k, :, :]
            b[l, :] .-= ratios .* b[k, :]
        end
    end
end

C = [1 2; 3 4]
C[:,1]



A = Float64.(cat([1 2; 3 4], [5 6; 7 8], [1 2; 3 4], dims=3))
b = Float64.([1 3 5; 2 4 6])

A,b

A[:, :, 1]


gem3d!(A, b)

A,b


# TODO: add and test
function tt(v, A)
    r = similar(A)
    @inbounds for j = 1:size(A,2)
        @simd for i = 1:size(A,1)
            r[i,j] = v[j] * A[i,j] # fixed a typo here!
        end
    end
    r
end

N = 3
coeffs = zeros(size(A, 3), size(A, 1))

for k in (N-1):-1:1
    coeffs[:,k] = (b[k,:])'
    for m in (k+1):(N-1)
        coeffs[:,k] = coeffs[:,k]-A[k,m,:] .* coeffs[:,m]
    end
    coeffs[:,k]=coeffs[:,k]./A[k,k,:];
end

coeffs


function back_substitution3d(A, b)
    row_count = size(A, 1)
    coeffs = zeros(size(A, 3), size(A, 1))
    for k in row_count:-1:1
        coeffs[:,k] = (b[k,:])'
        for m in (k+1):row_count
            coeffs[:,k] = coeffs[:,k]-A[k,m,:] .* coeffs[:,m]
        end
        coeffs[:,k]=coeffs[:,k]./A[k,k,:];
    end
    return coeffs
end


function multi_mat_div_gem(A, b)
    gem3d!(A, b)
    return back_substitution3d(A, b)
end


coeffs = multi_mat_div_gem(A, b)





coeffs[:, 1]


big = rand(1234, 1234)

##
using TimerOutputs
to = TimerOutput()

function foo()
   global to
   @timeit to "sleepytime" out = 1+1
   println(out)
   return nothing
end

foo()

print_timer(to)

@info "With debug timings disabled:"
TimerOutputs.disable_debug_timings(LAP_julia)
foo()

@info "With debug timings enabled:"
TimerOutputs.enable_debug_timings(LAP_julia)
foo()



## make spaghetti
using LAP_julia

img = gen_spaghetti((256,256), 0.15, 70, spread=20)
imgshow(img)


img = gen_one_spaghetti((200,200), 0.02, 10)
imgshow(img)

minimum(img)
maximum(img)

pl = Rotate(35);
img_processed = augment(img, pl);
imgshow(img_processed)


img_size = (200,200)
X = ones(img_size[1]) * collect(range(-1,1,length=img_size[2]))'
X_rot = augment(X, pl)


imgshow(X)
imgshow(X_rot)

img = rand(200,200)
pl = Either(1=>FlipX(), 1=>FlipY(), 2=>NoOp()) |>
            Rotate(0:360) |>
            ShearX(-5:5) * ShearY(-5:5) |>
            CropSize(165, 165) |>
            Zoom(1:0.05:1.2) |>
            Resize(64, 64)

img_new = augment(img, pl)
