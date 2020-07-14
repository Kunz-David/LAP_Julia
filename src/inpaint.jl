using ImageFiltering: kernelfactors, imfilter!, imfilter, centered

"""
    inpaint_nans!(flow::Flow)

Inpaint NaNs in `flow` using surrounding non-NaN values.

# References
Oliveira, Manuel & Bowen, Brian & McKenna, Richard & Chang, Yu-Sung. (2001). Fast Digital Image Inpainting.. 261-266.
[Here.](https://www.researchgate.net/publication/220903053_Fast_Digital_Image_Inpainting)

See also: [`showflow`](@ref), [`Flow`](@ref)
"""
function inpaint_nans!(flow::Flow)

    # make a mask where nans are 0 and the rest is 1
    mask = .!isnan.(real(flow))
    @views flow[.!mask] .= 0

    # init interpolation kernel
    ker = centered([0.5, 2, 0.5])
    kernf = kernelfactors((ker, ker))

    zero_in_mask = any(mask .== 0)

    # init tmp vars
    div_coef = similar(flow, Float16)
    estim_flow = similar(flow)

    # until there 0s left in mask
    while zero_in_mask

        imfilter!(div_coef, mask, kernf, "symmetric")
        imfilter!(estim_flow, mask .* flow, kernf, "symmetric")

        # identify the newly estimated
        newly_inpainted = ((div_coef .!= 0) .* (mask .== 0))

        # set newly inpainted
        @views flow[newly_inpainted] = estim_flow[newly_inpainted] ./ div_coef[newly_inpainted]

        # update mask
        @views mask[newly_inpainted] = ones(sum(newly_inpainted))

        # update while condition
        zero_in_mask = any(mask .== 0)
    end

    return nothing
end
