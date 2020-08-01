


## using Float64
timer=TimerOutput("sparse pflap psnr reg"); # fix
method_kwargs =Dict(:timer => timer, :display => false, :max_repeats => 3, :match_source_histogram => false)
flow_est, source_reg, timer, results = test_registration_alg(sparse_pflap_psnr, img, imgw, flow, method_kwargs=method_kwargs, timer=timer);










## using Float16
timer=TimerOutput("sparse pflap psnr float 16 reg"); # fix
method_kwargs =Dict(:timer => timer, :display => false, :max_repeats => 3, :match_source_histogram => false)
flow_est, source_reg, timer, results = test_registration_alg(LAP_julia.sparse_pflap_psnr_float16, img, imgw, flow, method_kwargs=method_kwargs, timer=timer);



filter_half_size = 15
filter_size = filter_half_size*2 +1
img_size = (257,257)

sigma = Float16.((filter_half_size + 2) / 4)
centered_inds = centered(-filter_half_size:filter_half_size)
gaus = KernelFactors.gaussian(sigma, filter_size)
kern16 = ImageFiltering.kernelfactors((gaus,gaus))


img16 = rand(Float16, img_size)
holder16 = rand(Float16, img_size)

@code_warntype imfilter!(holder16, img16, kern16)

@benchmark imfilter!($holder16, $img16, $kern16)

# --------->> for some reason 16 is slow to filter.


##

sigma = Float64.((filter_half_size + 2) / 4)
centered_inds = centered(-filter_half_size:filter_half_size)
gaus = KernelFactors.gaussian(sigma, filter_size)
kern64 = ImageFiltering.kernelfactors((gaus,gaus))

img64 = rand(Float64, img_size)
holder64 = rand(Float64, img_size)

@code_warntype imfilter!(holder64, img64, kern64)

@benchmark imfilter!($holder64, $img64, $kern64)



using LinearAlgebra

A = [1 1; 2 2]

b = [1, 1]

A\b
qr(A) \ b
qr(A, Val(true)) \ b
lu(A) \ b
svd(A) \ b


A = Float16[0.0 0.0; 0.0 0.0]
b = Float16[-0.009766, 0.0]
