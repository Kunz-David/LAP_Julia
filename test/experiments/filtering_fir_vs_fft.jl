using ImageFiltering
using LAP_julia: prepare_gaussian_filters


filter_speed_df = DataFrame(
    index = Int[],
    filter_alg = Symbol[],
    img_size = Int[],
    khs = Int[],
    bench = BenchmarkTools.Trial[])

kernel_half_sizes = cat(collect(1:51), collect(61:10:101), collect(141:40:501), collect(501:98:795), dims=1);
img_sizes = [100, 200, 400, 800, 1600, 3200];
filter_algs = [Algorithm.FIR, Algorithm.FFT]


##
df = deepcopy(filter_speed_df)
let index = 0
    for img_size in img_sizes

        rand_img = rand(img_size, img_size)

        # limit the max kernel size to 1/4 the image
        kernel_limit = img_size/4
        khs_modified = filter(x -> x <= kernel_limit, kernel_half_sizes)

        for khs in khs_modified

            forward_ker, backward_ker = prepare_gaussian_filters(khs)

            for alg in filter_algs
                filt_result = similar([1.0], (img_size, img_size))

                # do the filtering, timed
                bench = @benchmark imfilter!($alg, $filt_result, $rand_img, $forward_ker[1], "symmetric")

                index = index + 1
                println("at index: ", index, " img_size:", img_size, " khs: ", khs)
                push!(df, Dict(:index => index,
                               :alg => Symbol(alg),
                               :img_size => img_size,
                               :khs => khs,
                               :bench => bench))
            end
        end
    end
end


filt_result = similar([1.0], (img_size, img_size))
img = gen_chess(25,2)
alg = Algorithm.FIR
res = AbstractResource{alg}
forward_ker, backward_ker = prepare_gaussian_filters(13)
# do the filtering, timed
imfilter(Float64, rand_img, forward_ker[1], "symmetric", Algorithm.FIR)
bench = @benchmark imfilter!($res, $filt_result, $rand_img, $forward_ker[1], "symmetric")
