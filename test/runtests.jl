

using .LAP_julia
using Test
using FileIO, Images, Colors
using CSV

base_path = "/Users/MrTrololord/Google_Drive/cvut/bakalarka/anhir/"
dataset_path = "/Users/MrTrololord/Google_Drive/cvut/bakalarka/anhir/dataset/"
landmark_path = "/Users/MrTrololord/Google_Drive/cvut/bakalarka/anhir/landmarks/"
loc_table = CSV.read(base_path * "location_table.csv")

train_rows = loc_table[loc_table[:status] .== "training", :]

@testset "LAP_julia.jl" begin
    TEST_SIZE = 2
    @testset "Tests on random images $k" for k = 1:TEST_SIZE
        random_rows = train_rows[rand(1:size(train_rows, 1), TEST_SIZE), :]
        source_loc = (random_rows[k, :])[Symbol("Source image")]
        target_loc = (random_rows[k, :])[Symbol("Target image")]
        source_in = load(dataset_path * source_loc)
        target_in = load(dataset_path * target_loc)
        gray_source = Gray.(source_in)
        gray_target = Gray.(target_in)

        @testset "rescale_intensities" begin
            f((x,y)) = [x[:]; y[:]]
            @test maximum(f(LAP_julia.rescale_intensities(gray_source, gray_target))) == 255
            @test minimum(f(LAP_julia.rescale_intensities(gray_source, gray_target))) == 0
        end
    end
end


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
