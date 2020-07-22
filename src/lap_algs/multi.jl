
using ImageFiltering: Fill, KernelFactors.gaussian, centered, kernelfactors, imfilter!, padarray, Pad
using LinearAlgebra: qr
using PyPlot: Figure
using PyPlot
using BenchmarkTools

#TODO: add see section about timing in the docs to the api
"""
    polyfilter_lap(target::Image, source::Image; filter_num::Integer=3, max_repeats::Integer=1, display::Bool=true)

Find a transformation flow (complex displacment field), that transforms image `source` to image `target`.

# Arguments
- `target::Image`: the image we want `source` to look like.
- `source::Image`: warped image we want to transform into `target`.
- `filter_num::Integer=3`: the number of basis filters used in `single_lap` calls (so far only =3 implemented).
- `max_repeats::Integer=1`: the maximum number of times an iteration of one filter size can be repeated.
- `display::Bool=true`: use verbose prints and return an array of figures.
- `timer::TimerOutput=TimerOutput("polyfilter_lap")`: provide a timer which times certain blocks in the function.

# Outputs
- `flow::Flow`: is the complex vector field that transforms `source` closer to `target`.
- `source_reg::Image`: is the image `source` transformed by `flow`.
- [`figs::Matrix{Figure}`: is a 2D array of `PyPlot` Figures which shows the work of the algorithm at each iteration.
For each iteration there are 3 Figures in this order: 1) current `u_est`, 2) newest addition to `u_est` `Δ_u`, 3) current `source_reg`.]

# Describtion
Implements the basic concept of `Algorithm 2` from the [paper](http://www.ee.cuhk.edu.hk/~tblu/monsite/pdfs/gilliam1701.pdf) without some features.
It uses `single_lap` iteratively; in each iteration using the transformation estimated by `single_lap` to warp the `source` image closer to
the `target` image and then using this warped closer image as the source image in the next iteration, while using progressively smaller
`filter_half_sizes` to estimate even small and faster varying displacements.

See also: [`single_lap`](@ref), [`imgshow`](@ref),[`imgshowflow`](@ref), [`warp_imgshowflow`](@ref), [`Flow`](@ref)
"""
function polyfilter_lap(target::Image,
                        source::Image;
                        filter_num::Integer=3,
                        max_repeats::Integer=1,
                        display::Bool=true,
                        timer::TimerOutput=TimerOutput("polyfilter_lap"))

    @timeit_debug timer "setup" begin
        # convert images to floats
        source = Float64.(source)
        target = Float64.(target)

        # rescale images to have the whole [0, 1] spectrum.
        target, source = rescale_intensities(target, source)
        #NOTE: a histogram match might be a good idea (https://juliaimages.org/stable/function_reference/#Images.histmatch)

        # pad with zeros if sizes difer.
        target, source = pad_images(target, source)

        image_size = size(target)

        # set number of layers in the filter pyramid.
        level_count = floor(Int64, log2(minimum(size(target))/8)+1)+1

        @assert ((2^(level_count)+1) <= minimum(image_size)) "level number results in a filter larger than the size of the input images."

        # displacement init.
        u_est = zeros(ComplexF64, size(target))

        # filter half sizes array eg. [16, 8, 4, 2, 1]
        half_size_pyramid = Int64.(2 .^ range(level_count-1, stop=0, length=level_count))


        if display
            num_plots = 4
            figs = Array{Figure}(undef, level_count, max_repeats*num_plots)
        end

        source_reg = source
    end

    for level in 1:level_count

        @timeit_debug timer "single filter pyramid level" begin
            if display
                println("###################")
                println("ITERATION: ", level)
                println("filter_half_size: ", half_size_pyramid[level])
            end

            filter_half_size = half_size_pyramid[level]::Int
            window_size = UInt8.(2 .* filter_half_size .* (1, 1) .+ 1)
            window_half_size = Int64.((window_size .- 1) ./ 2)

            for iter_repeat in 1:max_repeats
                @timeit_debug timer "LAP" begin
                    Δ_u = single_lap(target, source_reg, filter_half_size, window_size, filter_num, timer=timer)
                end

                # USE INPAINTING TO CORRECT U_EST:
                # 1) replicate borders
                @timeit_debug timer "inpainting" begin
                    @timeit_debug timer "replicating borders" begin
                        middle_vals = Δ_u[window_half_size[1]+1:end-window_half_size[1],
                                          window_half_size[2]+1:end-window_half_size[2]]
                        Δ_u = parent(padarray(middle_vals, Pad(:replicate, window_half_size...)))
                    end # "replication borders"
                        # 2) inpaint middle missing values
                    if all(isnan.(real(Δ_u))) # if all are NaNs
                        Δ_u = zeros(ComplexF64, image_size)
                    elseif any(isnan.(Δ_u)) # if some are NaNs
                        inpaint_nans!(Δ_u)
                    end
                end # "inpainting"

                # SMOOTH U_EST WITH A GAUSSIAN FILTER:
                @timeit_debug timer "smoothing" begin
                    Δ_u = smooth_with_gaussian!(Δ_u, window_half_size)
                end

                # add the Δ_u to my u_est
                u_est = u_est + Δ_u

                @timeit_debug timer "image interpolation" begin
                    # linear interpolation
                    source_reg = warp_img(source, real(u_est), imag(u_est))
                end # "interpolation"

                if display
                    # @assert eltype(Δ_u) <: Complex
                    figs[level, 1] = showflow(u_est, figtitle="U_EST (Level: " * string(level) * "/" * string(level_count) * ")")
                    figs[level, 2] = showflow(Δ_u, figtitle="Δ_U (Level: " * string(level) * "/" * string(level_count) * ")")
                    figs[level, 3] = imgshow(source_reg, figtitle="SOURCE_REG (Level: " * string(level) * "/" * string(level_count) * ")")
                end
            end
        end # "single filter pyramid level"
    end

    # if display; return u_est, source_reg, figs; end
    return u_est, source_reg
end


function sparse_pflap(target::Image,
                      source::Image;
                      filter_num::Integer=3,
                      max_repeats::Integer=1,
                      display=true,
                      point_count::Int=35,
                      spacing::Int=35,
                      timer::TimerOutput=TimerOutput("sparse_pflap"))

    @timeit_debug timer "setup" begin
        # convert images to floats
        source = Float64.(source)
        target = Float64.(target)


        #NOTE: a histogram match might be a good idea (https://juliaimages.org/stable/function_reference/#Images.histmatch)
        # adjust the source histogram to be like target
        @timeit_debug timer "hist match" begin
            source = adjust_histogram(source, Matching(targetimg = target)) # takes 8ms on average
        end

        # pad with zeros if sizes difer.
        target, source = pad_images(target, source)

        image_size = size(target)

        # set number of layers in the filter pyramid.
        level_count = floor(Int64, log2(minimum(size(target))/8)+1)+1

        @assert ((2^(level_count)+1) <= minimum(image_size)) "level number results in a filter larger than the size of the input images."

        # displacement init.
        u_est = Array{Complex{Float64},2}(undef, image_size...)

        # filter half sizes array eg. [16, 8, 4, 2, 1]
        half_size_pyramid = Int64.(2 .^ range(level_count-1, stop=0, length=level_count))

        if display
            num_plots = 5
            figs = Array{Figure}(undef, level_count, max_repeats*num_plots)
        end

        source_reg = source

        #TODO edit
        # get edge points
        mask = falses(size(target))
        fhs = 3
        mask[fhs+1:end-fhs, fhs+1:end-fhs] .= true
        @timeit_debug timer "find edge points" begin
            inds = find_edge_points(target, spacing=spacing, number=point_count, mask=mask)
        end
        if display
            println("ind count: ", length(inds))
        end

    end # "setup"

    for level in 1:level_count

        @timeit_debug timer "single filter pyramid level" begin
            if display
                println("###################")
                println("ITERATION: ", level)
                println("filter_half_size: ", half_size_pyramid[level])
            end

            fhs = half_size_pyramid[level]
            window_size = Int64.(2 .* fhs .+ (1, 1))

            if display
                println("inds: ", length(inds))
            end

            for iter_repeat in 1:max_repeats
                @timeit_debug timer "single lap at points" begin
                    new_estim_at_inds = single_lap_at_points(target, source_reg, fhs, window_size, inds, timer=timer, display=display)
                end

                # interpolate flow
                if all(isnan, new_estim_at_inds)
                    @timeit_debug timer "interpolate flow" begin
                        if display
                            println("\tIMPORTANT: all new estim vectors are NaN.")
                        end
                    end
                else
                    @timeit_debug timer "interpolate flow" begin
                        # extract u_est_vecs from the previous estimate and add them to the newly estimated ones
                        # then interpolate a new flow.
                        if display
                            println("non NaN new estim flow vector count: ", count(!isnan, new_estim_at_inds))
                        end
                        u_est_at_inds = map(ind -> u_est[ind], inds) .+ new_estim_at_inds
                        u_est = interpolate_flow(u_est_at_inds, inds, size(u_est))
                    end
                end

                @assert any(.!isnan.(real(u_est)))

                # linear interpolation
                @timeit_debug timer "interpolate image" begin
                    source_reg = warp_img(source, real(u_est), imag(u_est))
                end

                if display
                    figs[level, 1] = showflow(u_est, figtitle="U_EST (Level: " * string(level) * "/" * string(level_count) * ")")
                    # figs[level, 2] = showflow(Δ_u, figtitle="Δ_U (Level: " * string(level) * "/" * string(level_count) * ")")
                    figs[level, 3] = imgshow(source_reg, figtitle="SOURCE_REG (Level: " * string(level) * "/" * string(level_count) * ")")
                    figs[level, 4] = showflow(create_sparse_flow_from_sparse(new_estim_at_inds, inds, size(u_est)), disp_type=:sparse, figtitle="SPARSE Δ_U (Level: " * string(level) * "/" * string(level_count) * ")")
                    figs[level, 5] = showflow(create_sparse_flow_from_full(u_est, inds), figtitle="u_est_at_points (Level: " * string(level) * "/" * string(level_count) * ")"); PyPlot.scatter([ind[2] for ind in inds], [ind[1] for ind in inds], marker = :x); gcf()
                end
            end
        end # "single filter pyramid level"
    end

    # here
    if display; return u_est, source_reg, figs; end
    return u_est, source_reg
end


function sparse_pflap_save(target::Image,
                                  source::Image;
                                  filter_num::Integer=3,
                                  max_repeats::Integer=1,
                                  display::Bool=true,
                                  point_count::Int=25,
                                  spacing::Int=40)

    # convert images to floats
    source = Float64.(source)
    target = Float64.(target)

    # rescale images to have the whole [0, 1] spectrum.
    target, source = rescale_intensities(target, source)
    #NOTE: a histogram match might be a good idea (https://juliaimages.org/stable/function_reference/#Images.histmatch)

    # pad with zeros if sizes difer.
    target, source = pad_images(target, source)

    image_size = size(target)

    # set number of layers in the filter pyramid.
    level_count = floor(Int64, log2(minimum(size(target))/8)+1)+1

    @assert ((2^(level_count)+1) <= minimum(image_size)) "level number results in a filter larger than the size of the input images."

    # displacement init.
    u_est = zeros(image_size) .+ zeros(image_size) .* im
    Δ_u = zeros(image_size) .+ zeros(image_size) .* im

    # filter half sizes array eg. [16, 8, 4, 2, 1]
    half_size_pyramid::Array{Int64,1} = 2 .^ range(level_count-1, stop=0, length=level_count)

    # at what filter size change the interpolation strategy
    interpol_change_index = findfirst(x -> x == 2, half_size_pyramid)

    if display
        num_plots = 5
        figs = Array{Figure}(undef, level_count, max_repeats*num_plots)
    end

    source_reg = source

    for level in 1:level_count

        if display
            println("###################")
            println("ITERATION: ", level)
            println("filter_half_size: ", half_size_pyramid[level])
        end

        fhs = half_size_pyramid[level]::Int
        window_size::Tuple{Int64, Int64} = 2 .* fhs .* (1, 1) .+ 1
        window_half_size::Tuple{Int64, Int64} = (window_size .- 1) ./ 2

        mask = parent(padarray(trues(size(target).-(2*fhs, 2*fhs)), Fill(false, (fhs, fhs), (fhs, fhs))))
        inds = find_edge_points(target, spacing=spacing, number=point_count, mask=mask)
        println("inds: ", length(inds))

        for iter_repeat in 1:max_repeats

            Δ_u_at_points = single_lap_at_points(target, source_reg, fhs, window_size, inds, filter_num)

            Δ_u_interpolated = interpolate_flow(Δ_u_at_points, inds)

            # USE INPAINTING TO CORRECT U_EST:
            # 1) replicate borders
            # middle_vals = Δ_u[window_half_size[1]+1:end-window_half_size[1],
            #                   window_half_size[2]+1:end-window_half_size[2]]
            # Δ_u = parent(padarray(middle_vals, Pad(:replicate, window_half_size...)))
            # # 2) inpaint middle missing values
            if all(isnan.(real(Δ_u_interpolated))) # if all are NaNs
                Δ_u = zeros(image_size) .+ zeros(image_size) .* im
            else
                Δ_u = Δ_u_interpolated
            end
            # # SMOOTH U_EST WITH A GAUSSIAN FILTER:
            # Δ_u = LAP_julia.smooth_with_gaussian(Δ_u, window_half_size)

            # add the Δ_u to my u_est
            u_est = u_est + Δ_u

            # linear interpolation
            source_reg = warp_img(source, real(u_est), imag(u_est))

            if display
                figs[level, 1] = showflow(u_est, figtitle="U_EST (Level: " * string(level) * "/" * string(level_count) * ")")
                figs[level, 2] = showflow(Δ_u, figtitle="Δ_U (Level: " * string(level) * "/" * string(level_count) * ")")
                figs[level, 3] = imgshow(source_reg, figtitle="SOURCE_REG (Level: " * string(level) * "/" * string(level_count) * ")")
                figs[level, 4] = showflow(Δ_u_at_points, disp_type=:sparse, figtitle="SPARSE Δ_U (Level: " * string(level) * "/" * string(level_count) * ")")
                figs[level, 5] = showflow(Δ_u, figtitle="Δ_U with points (Level: " * string(level) * "/" * string(level_count) * ")"); PyPlot.scatter([ind[2] for ind in inds], [ind[1] for ind in inds], marker = :x); gcf()
            end
        end

        if display
            println("###################")
        end
    end

    # here
    if display; return u_est, source_reg, figs, Δ_u; end
    return u_est, source_reg
end
