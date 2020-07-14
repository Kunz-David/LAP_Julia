
using Test
using FileIO, Images, Colors
using CSV


# !!! commented out because of local repository of test images. !!!

# base_path = "/Users/MrTrololord/Google_Drive/cvut/bakalarka/anhir/"
# dataset_path = "/Users/MrTrololord/Google_Drive/cvut/bakalarka/anhir/dataset/"
# landmark_path = "/Users/MrTrololord/Google_Drive/cvut/bakalarka/anhir/landmarks/"
# loc_table = CSV.read(base_path * "location_table.csv")
#
# train_rows = loc_table[loc_table[:status] .== "training", :]
#
# @testset "LAP_julia.jl" begin
#     TEST_SIZE = 2
#     @testset "Tests on random images $k" for k = 1:TEST_SIZE
#         random_rows = train_rows[rand(1:size(train_rows, 1), TEST_SIZE), :]
#         source_loc = (random_rows[k, :])[Symbol("Source image")]
#         target_loc = (random_rows[k, :])[Symbol("Target image")]
#         source_in = load(dataset_path * source_loc)
#         target_in = load(dataset_path * target_loc)
#         gray_source = Gray.(source_in)
#         gray_target = Gray.(target_in)
#
#         @testset "rescale_intensities" begin
#             both((x,y)) = [x[:]; y[:]]
#             @test maximum(both(LAP_julia.helpers.rescale_intensities(gray_source, gray_target))) == 1
#             @test minimum(both(LAP_julia.helpers.rescale_intensities(gray_source, gray_target))) == 0
#         end
#     end
# end

@testset "Registration Algorithms" begin
    TimerOutputs.enable_debug_timings(LAP_julia.registration)
    img, imgw, flow = gen_init(:lena, flow_args=[15, 140])

    @testset "PFLAP" begin
        timer = TimerOutput("Registration");
        @timeit timer "polyfilter lap" begin
            flow_est, source_reg = polyfilter_lap(img, imgw, display=false, timer=timer)
        end
        print_timer(timer)

        flow_est, source_reg = polyfilter_lap(img, imgw, display=false)
        showflow(flow)
        showflow(flow_est)
        @test angle_mae(flow, flow_est) < 10










@testset "helpers" begin
    @testset "pad_images $k, $l, $m, $n" for k = 1:3, l = 1:3, m = 1:3, n = 1:3
        im1 = ones(k, l)
        im2 = ones(m, n)
        @test size(LAP_julia.pad_images(im1, im2)[1]) <= (k+m-1, l+n-1)
        @test size(LAP_julia.pad_images(im1, im2)[1]) == ((k<m) ? m : k, (l<n) ? n : l)
    end
end

@testset "inpaint" begin
    u = ones(10,10) .+ 2im .* ones(10,10);
    nan_u = u
    nan_u[3:8, 3:8] .= NaN .+ NaN .* 1im;
    LAP_julia.inpaint.inpaint_nans!(nan_u)
    @test all(nan_u .== u)
    @test any(isnan.(nan_u)) == false
end
