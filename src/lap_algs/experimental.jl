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
function sparse_pflap_psnr(target::Image,
                           source::Image;
                           filter_count::Integer=3,
                           max_repeats::Integer=3,
                           display::Bool=false,
                           point_count::Int=500,
                           spacing::Int=13,
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
        u_est = zeros(ComplexF64, size(target))

        # filter half sizes array eg. [16, 8, 4, 2, 1]
        half_size_pyramid = Int64.(2 .^ range(level_count-1, stop=0, length=level_count))

        if display
            num_plots = 4
            figs = Array{Figure}(undef, level_count, max_repeats, num_plots)
        end

        psnr_threshold = 0.3
        psnr_before = 0
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
                    if display
                        println("\tnon border inds: $(length(current_inds))/$(length(inds)) ")
                    end
                    if level != 1
                        current_inds = LAP_julia.mapped_out(current_inds, target)
                    end
                end
                if display
                    println("inner loop iter: $iter_repeat")
                    println("\tnon mapped out inds: $(length(current_inds))/$(length(inds)) ")
                end
                @timeit_debug timer "single lap at points" begin
                    new_estim_at_inds, current_inds = single_lap_at_points(target, source_reg, fhs, window_size, current_inds, timer=timer, display=display)
                end

                # interpolate flow
                if all(isnan, new_estim_at_inds)
                    @timeit_debug timer "interpolate flow" begin
                        if display
                            println("\tIMPORTANT: all new estim vectors are NaN.")
                        end
                    end
                    break
                else
                    @timeit_debug timer "flow interpolation" begin
                        # extract u_est_vecs from the previous estimate and add them to the newly estimated ones
                        # then interpolate a new flow.
                        if display
                            println("\tnon NaN new estim flow vector count: ", length(new_estim_at_inds))
                        end
                        u_est_at_inds = map(ind -> u_est[ind], current_inds) .+ new_estim_at_inds
                        # @assert any(.!isnan.(u_est_at_inds)) "new estim nans", count(isnan, u_est_at_inds), length(u_est_at_inds)
                        u_est_adept = interpolate_flow(u_est_at_inds, current_inds, size(u_est))
                    end
                end


                # if the interpolation failes and returns u_est_adept full of NaNs go next level.
                if any(isnan.(u_est_adept))
                    if display
                        println("\tIMPORTANT: u_est_adept full of NaNs, skipping")
                    end
                    break
                end

                # if the interpolation estimated a flow 20% larger than the fhs then go to next level
                # if any(vec_len.(u_est_adept) .> (fhs*1.3))
                #     if display
                #         println("\tIMPORTANT: u_est_adept is above size threshold, skipping (max: $(maximum(vec_len.(u_est_adept))), fhs: $fhs)")
                #     end
                #     break
                # end
                @assert any(.!isnan.(u_est_adept)) count(isnan, u_est_adept), length(u_est_adept)

                # linear interpolation
                @timeit_debug timer "image interpolation" begin
                    # source_reg = warp_img(source, -real(u_est_adept), -imag(u_est_adept), target)
                    source_reg = warp_img(source, -real(u_est_adept), -imag(u_est_adept))
                end

                # calculate the psnr after warping the source by u_est_adept
                psnr_after = assess_psnr(source_reg, target)

                # compare psnr before and after to see if there was an improvement larger than psnr_threshold
                if psnr_after - psnr_before > psnr_threshold
                    # nice, the improvement is large -> try to run again with the same filter size.
                    u_est = u_est_adept
                    psnr_before = psnr_after
                    if display
                        println("\tstatus: success, repeating filter size")
                    end
                else
                    # the improvement isn't larger than psnr_threshold -> next filter size
                    # check to see if it was even an improvement:
                    if psnr_after - psnr_before > 0
                        # it was an improvement, save the variables:
                        u_est = u_est_adept
                        psnr_before = psnr_after
                        if display
                            println("\tstatus: success, next filter size")
                            add_figs_sparse_pflap(figs, level, iter_repeat, u_est, new_estim_at_inds, source_reg, level_count, max_repeats, current_inds, fhs)
                        end
                    else
                        # it wasn't an improvement, skip
                        if display
                            println("\tstatus: failure, trashing")
                        end
                    end
                    break # go for next filter size.
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
