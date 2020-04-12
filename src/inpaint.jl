module inpaint

export inpaint_nans!

using ImageFiltering: kernelfactors, imfilter!, imfilter, centered

"""
    function inpaint_nans(u)

# input:
- `u` ... complex matrix containing some NaNs

# output:
- `u_out` ... complex matrix *not* containing NaNs

# source:
- (https://www.researchgate.net/publication/220903053_Fast_Digital_Image_Inpainting)
"""
function inpaint_nans!(u)

    # make a mask where nans are 0 and the rest is 1
    mask = .!isnan.(real(u))
    @views u[.!mask] .= 0

    # init interpolation kernel
    ker = centered([0.5, 2, 0.5])
    kernf = kernelfactors((ker, ker))

    zero_in_mask = any(mask .== 0)

    # init tmp vars
    div_coef = similar(u, Float16)
    estim_flow = similar(u)

    # until there 0s left in mask
    while zero_in_mask

        imfilter!(div_coef, mask, kernf, "symmetric")
        imfilter!(estim_flow, mask .* u, kernf, "symmetric")

        # identify the newly estimated
        newly_inpainted = ((div_coef .!= 0) .* (mask .== 0))

        # set newly inpainted
        @views u[newly_inpainted] = estim_flow[newly_inpainted] ./ div_coef[newly_inpainted]

        # update mask
        @views mask[newly_inpainted] = ones(sum(newly_inpainted))

        # update while condition
        zero_in_mask = any(mask .== 0)
    end

    return u
end

end #module
