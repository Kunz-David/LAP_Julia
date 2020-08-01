
using ImageFiltering: Fill, KernelFactors.gaussian, centered, kernelfactors, imfilter!, padarray, Pad
using LinearAlgebra: qr
using PyPlot: Figure
using PyPlot, BenchmarkTools
using Images: assess_psnr

#TODO: add see section about timing in the docs to the api
"""
    pflap(args; kwargs)

Find a transformation flow (complex displacment field), that transforms image `source` to image `target`.
Returns the transformation a [`Flow`](@ref) and the registered `source` image.

# Arguments
- `target::Image`: the image we want `source` to look like.
- `source::Image`: warped image we want to transform into `target`.

# Keyword Arguments
- `filter_count::Integer=3`: the number of basis filters used in `single_lap` calls (so far only =3 implemented).
- `max_repeats::Integer=1`: the maximum number of times an iteration of one filter size can be repeated.
- `display::Bool=true`: use verbose prints and return an array of figures.
- `timer::TimerOutput=TimerOutput("pflap")`: provide a timer which times certain blocks in the function.
- `match_source_histogram::Bool=true`: choose to match the `source` image histogram to the `target` image histogram.
- `rescale_intensities::Bool=false`: choose to rescale! intensities of the images to the whole `[0, 1]` spectrum.

# Outputs
- `flow::Flow`: is the complex vector field that transforms `source` closer to `target`.
- `source_reg::Image`: is the image `source` transformed by `flow`.
- [`figs::Matrix{Figure}`: is a 2D array of `PyPlot` Figures which shows the work of the algorithm at each iteration. For each iteration there are 3 Figures in this order: 1) current `u_est`, 2) newest addition to `u_est` `Δ_u`, 3) current `source_reg`.]

# Describtion
Implements the basic concept of `Algorithm 2` from the [paper](http://www.ee.cuhk.edu.hk/~tblu/monsite/pdfs/gilliam1701.pdf) without some features.
It uses `single_lap` iteratively; in each iteration using the transformation estimated by `single_lap` to warp the `source` image closer to
the `target` image and then using this warped closer image as the source image in the next iteration, while using progressively smaller
`filter_half_sizes` to estimate even small and faster varying displacements.

## Note:
If `display=true` is chosen the function will also return an array of figures demonstrating the
work of the registration.

See also: [`single_lap`](@ref), [`imgshow`](@ref),[`imgshowflow`](@ref), [`warp_imgshowflow`](@ref), [`Flow`](@ref)
"""
function pflap(target::Image,
                        source::Image;
                        filter_count::Integer=3,
                        max_repeats::Integer=3,
                        display::Bool=true,
                        timer::TimerOutput=TimerOutput("pflap"),
                        prefilter::Bool=false,
                        match_source_histogram::Bool=false,
                        rescale_intensities::Bool=false)

    @timeit_debug timer "setup" begin
        # convert images to floats
        source = Float64.(source)
        target = Float64.(target)

        # rescale images to have the whole [0, 1] spectrum.
        if rescale_intensities
            target, source = rescale!(target, source)
        end

        # match histogram
        if match_source_histogram
            @timeit_debug timer "hist match" begin
                source = adjust_histogram(source, Matching(targetimg = target)) # takes 8ms on average
            end
        end

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

        # threshold above which we repeat the iteration.
        psnr_threshold=0.3

        if display
            num_plots = 4
            figs = Array{Figure}(undef, level_count, max_repeats, num_plots)
        end

        psnr_before = 0

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
            window_size = Int64.(2 .* filter_half_size .* (1, 1) .+ 1)
            window_half_size = Int64.((window_size .- 1) ./ 2)

            if prefilter
                @timeit_debug timer "prefiltering" begin
                    target_inner = highpass_image(target, filter_half_size)
                end
            else
                target_inner = target
            end

            for iter_repeat in 1:max_repeats

                if display
                    println("inner loop iter: $iter_repeat")
                end

                if prefilter
                    @timeit_debug timer "prefiltering" begin
                        source_reg_inner = highpass_image(source_reg, filter_half_size)
                    end
                else
                    source_reg_inner = source_reg
                end

                @timeit_debug timer "LAP" begin
                    Δ_u = single_lap(target_inner, source_reg_inner, filter_half_size, window_size, filter_count, timer=timer)
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
                        if display
                            println("NaN count: $(count(isnan, Δ_u))")
                        end
                        inpaint_nans!(Δ_u)
                    end
                end # "inpainting"

                # SMOOTH U_EST WITH A GAUSSIAN FILTER:
                @timeit_debug timer "smoothing" begin
                    Δ_u = smooth_with_gaussian!(Δ_u, window_half_size)
                end

                # add the Δ_u to my u_est
                u_est_adept = u_est + Δ_u

                @timeit_debug timer "image interpolation" begin
                    # linear interpolation
                    source_reg = warp_img(source, -real(u_est_adept), -imag(u_est_adept))
                    # source_reg = warp_img(source, -real(u_est_adept), -imag(u_est_adept), target, border_strat=:fill_img)
                end # "interpolation"

                if prefilter
                    @timeit_debug timer "prefiltering" begin
                        source_reg_inner = highpass_image(source_reg, filter_half_size)
                    end
                else
                    source_reg_inner = source_reg
                end

                # calculate the psnr after warping the source by u_est_adept
                psnr_after = assess_psnr(source_reg_inner, target_inner)

                # compare psnr before and after to see if there was an improvement larger than psnr_threshold
                if psnr_after - psnr_before > psnr_threshold
                    # nice, the improvement is large -> try to run again with the same filter size.
                    u_est = u_est_adept
                    psnr_before = psnr_after
                    if display
                        println("\t status: success, repeating filter size")
                    end
                else
                    # the improvement isn't larger than psnr_threshold -> next filter size
                    # check to see if it was even an improvement:
                    if psnr_after - psnr_before > 0
                        # it was an improvement, save the variables:
                        u_est = u_est_adept
                        psnr_before = psnr_after
                        if display
                            println("\t status: success, next filter size")
                            add_figs_pflap(figs, level, iter_repeat, u_est, Δ_u, source_reg, level_count, max_repeats, filter_half_size)
                        end
                    else
                        # it wasn't an improvement, skip
                        if display
                            println("\t status: failure, trashing")
                        end
                    end
                    break # go for next filter size.
                end

                if display
                    add_figs_pflap(figs, level, iter_repeat, u_est, Δ_u, source_reg, level_count, max_repeats, filter_half_size)
                end
            end
        end # "single filter pyramid level"
    end

    if display; return u_est, source_reg, figs; end
    return u_est, source_reg
end


"""
    add_figs_pflap(figs, level, iter_repeat, u_est, Δ_u, source_reg, level_count, max_repeats)

Fill the figs array specifically for the `pflap` method.
"""
function add_figs_pflap(figs, level, iter_repeat, u_est, Δ_u, source_reg, level_count, max_repeats, fhs)
    title_string = "(Level: $(string(level))/$(string(level_count))), (Iter: $(string(iter_repeat))/$(string(max_repeats))), (fhs: $fhs)"
    figs[level, iter_repeat, 1] = showflow(u_est, figtitle="U_EST " * title_string)
    figs[level, iter_repeat, 2] = showflow(Δ_u, figtitle="Δ_U " * title_string)
    figs[level, iter_repeat, 3] = imgshow(source_reg, figtitle="SOURCE_REG " * title_string)
end

"""
    sparse_pflap(args; kwargs)

Find a transformation flow (complex displacment field), that transforms image `source` to closer to image `target`.
Returns the transformation a [`Flow`](@ref) and the registered `source` image.

# Arguments:
- `target::Image`: target/fixed grayscale image.
- `source::Image`: source/moving grayscale image.

# Keyword Arguments:
- `filter_count::Integer=3`: the number of basis filters used (so far only =3 implemented).
- `timer::TimerOutput=TimerOutput("sparse pflap")`: provide a timer which times certain blocks in the function.
- `display::Bool=false`: verbose and debug prints.
- `point_count::Int=500`: the number of points attempted to be found at the edges of the target image.
- `spacing::Int=10`: the minimal distance between two points.
- `match_source_histogram::Bool=true`: choose to match the `source` image histogram to the `target` image histogram.

## Note:
If `display=true` is chosen the function will also return an array of figures demonstrating the
work of the registration.

See also: [`sparse_lap`](@ref), [`pflap`](@ref), [`showflow`](@ref), [`single_lap`](@ref), [`sparse_pflap`](@ref).
"""
function sparse_pflap(target::Image,
                      source::Image;
                      filter_count::Integer=3,
                      max_repeats::Integer=1,
                      display::Bool=false,
                      point_count::Int=500,
                      spacing::Int=10,
                      timer::TimerOutput=TimerOutput("sparse pflap"),
                      match_source_histogram::Bool=true,
                      rescale_intensities::Bool=false)

    @timeit_debug timer "setup" begin
        # convert images to floats
        source = Float64.(source)
        target = Float64.(target)

        # rescale images to have the whole [0, 1] spectrum.
        if rescale_intensities
            target, source = rescale!(target, source)
        end

        # adjust the source histogram to be like target
        if match_source_histogram
            @timeit_debug timer "hist match" begin
                source = adjust_histogram(source, Matching(targetimg = target)) # takes 8ms on average
            end
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
            num_plots = 4
            figs = Array{Figure}(undef, level_count, max_repeats, num_plots)
        end

        source_reg = source

        #TODO edit
        # get edge points
        mask = falses(size(target))
        fhs = 3
        mask[fhs+1:end-fhs, fhs+1:end-fhs] .= true
        @timeit_debug timer "find edge points" begin
            inds = find_edge_points(target, spacing=spacing, point_count=point_count, mask=mask)
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

            for iter_repeat in 1:max_repeats

                @timeit_debug timer "filter inds" begin
                    current_inds = filter_border_inds(inds, size(target), fhs)
                end
                if display
                    println("inner loop iter: $iter_repeat")
                    println("\tcurrent inds: ", length(current_inds))
                end
                @timeit_debug timer "single lap at points" begin
                    new_estim_at_inds, current_inds = single_lap_at_points(target, source_reg, fhs, window_size, current_inds, timer=timer, display=display)
                end

                # interpolate flow
                if all(isnan, new_estim_at_inds)
                    @timeit_debug timer "flow interpolation" begin
                        if display
                            println("\tIMPORTANT: all new estim vectors are NaN.")
                        end
                    end
                else
                    @timeit_debug timer "flow interpolation" begin
                        # extract u_est_vecs from the previous estimate and add them to the newly estimated ones
                        # then interpolate a new flow.
                        if display
                            println("\tnon NaN new estim flow vector count: ", length(new_estim_at_inds))
                        end
                        u_est_at_inds = map(ind -> u_est[ind], current_inds) .+ new_estim_at_inds
                        @assert any(.!isnan.(u_est_at_inds)) "new estim nans", count(isnan, u_est_at_inds), length(u_est_at_inds)
                        u_est = interpolate_flow(u_est_at_inds, current_inds, size(u_est))
                    end
                end

                if any(isnan.(u_est))
                    if display
                        println("\tIMPORTNT: interpolation is full of NaNs, skip")
                    end
                    continue;
                end
                # @assert any(.!isnan.(u_est)) count(isnan, u_est), length(u_est)

                # linear interpolation
                @timeit_debug timer "image interpolation" begin
                    # source_reg = warp_img(source, -real(u_est), -imag(u_est))
                    source_reg = warp_img(source, -real(u_est), -imag(u_est), target)
                end

                if display
                    add_figs_sparse_pflap(figs, level, iter_repeat, u_est, new_estim_at_inds, source_reg, level_count, max_repeats, current_inds, fhs)
                end
            end
        end # "single filter pyramid level"
    end

    if display; return u_est, source_reg, figs; end
    return u_est, source_reg
end

function add_figs_sparse_pflap(figs, level, iter_repeat, u_est, new_estim_at_inds, source_reg, level_count, max_repeats, current_inds, fhs)
    level_iter_fhs = "(Level: $(string(level))/$(string(level_count))), (Iter: $(string(iter_repeat))/$(string(max_repeats))), (fhs: $fhs)"
    figs[level, iter_repeat, 1] = showflow(u_est, figtitle="U_EST " * level_iter_fhs)
    figs[level, iter_repeat, 2] = imgshow(source_reg, figtitle="SOURCE_REG " * level_iter_fhs, origin_left_bot=true)
    figs[level, iter_repeat, 3] = showflow(create_sparse_flow_from_sparse(new_estim_at_inds, current_inds, size(u_est)), disp_type=:sparse, figtitle="SPARSE Δ_U " * level_iter_fhs)
    figs[level, iter_repeat, 4] = showflow(create_sparse_flow_from_full(u_est, current_inds), figtitle="u_est_at_points " * level_iter_fhs)
end
