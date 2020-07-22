using LAP_julia: inds_to_points, asses_flow_quality

## testing

flow_size = size(flow)

## method
img, imgw, flow = gen_init()
timer=TimerOutput("reg alg: sparse lap");
method_kwargs = Dict(:timer => timer, :point_count => 200, :spacing => 15)
flow_est, source_reg, timer, results, (estim_at_inds, inds) = test_registration_alg(sparse_lap, img, imgw, flow, [25, (51,51)], method_kwargs, timer=timer)

showflow(flow)
showflow(flow_est)
##

estim_flow = interpolate_flow_quad(estim_at_inds, inds, size(flow))
estim_flow_rbf = interpolate_flow(estim_at_inds, inds, size(flow))

asses_flow_quality(flow, estim_flow)


asses_flow_quality(flow, estim_flow_rbf)

showflow(flow, figtitle="original")
showflow(estim_flow, figtitle="quad estimate")
showflow(estim_flow_rbf, figtitle="rbf estimate")


estim_flow = interpolate_flow_quad(estim_at_inds, inds, size(flow))
estim_flow_rbf = interpolate_flow(estim_at_inds, inds, size(flow))

asses_flow_quality(flow, estim_flow)
asses_flow_quality(flow, estim_flow_rbf)


showflow(estim_flow)
showflow(flow)
showflow(flow_truth .- estim_flow)


showflow(rbf_flow)
asses_flow_quality(flow, rbf_flow)
showflow(flow_truth .- rbf_flow)


@benchmark interpolate_flow_quad(estim_at_inds, inds, size(flow))


## benchmark both

rand_inds = rand(CartesianIndices(size(flow)), 15)
lsdkf = interpolate_flow_quad(map(x -> flow[x], rand_inds), rand_inds, size(flow))
flow
showflow(lsdkf)
showflow(flow)



# function make_flow_at_ind(ind, a, b, c, d, e, f)
#     return a +
#            b * ind[1] +
#            c * ind[2] +
#            d * ind[1]^2 +
#            e * ind[2]^2 +
#            f * ind[1]*ind[2]
# end


# function make_flow(coeffs, flow_size)
#     return reduce((coeff, first) -> coeff .* ones(flow_size) .+ first, coeffs; init=0)
# end
