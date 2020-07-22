using TimerOutputs


"""
    test_registration_alg(args, kwargs)

Test a registration function `reg_fun` by timing it and comparing its outputs to the ground truth `flow` and the `target` image.

# Arguments:
- `reg_fun`: registration function.
- `target::Image`: target image.
- `source::Image`: source image to be warped closer to `target`.
- `flow::Flow`: truth displacement flow warping `target` to `source`.
- `reg_fun_args=[]`: arguments for the registration function.
- `reg_fun_kwargs::Dict=Dict()`: keyword arguments for the registration function.

# Keyword Arguments:
- `display::Bool=true`: print debugging info and test results.
- `timer::TimerOutput=TimerOutput("blank")`: timer to for certain blocks of code in the `reg_fun`.

# Outputs:
- `flow_est`: estimated transformation flow.
- `source_reg`: `source` image registered.
- `timer`: timings of the `reg_fun` insides.
- `results`: Dict with results of the quality of the registration.
- [`(outputs)`]: other outputs of the `reg_fun` besides `flow_est` and `source_reg`.

See also: [`lap`](@ref), [`sparse_lap`](@ref),[`polyfilter_lap`](@ref), [`sparse_pflap`](@ref), [`Flow`](@ref)
"""
function test_registration_alg(reg_fun,
                               target::Image,
                               source::Image,
                               flow::Flow,
                               reg_fun_args=[],
                               reg_fun_kwargs::Dict=Dict();
                               display::Bool=true,
                               timer::TimerOutput=TimerOutput("blank"))
    TimerOutputs.enable_debug_timings(LAP_julia)

    @timeit timer timer.name begin
        outputs = reg_fun(target, source, reg_fun_args...; timer=timer, reg_fun_kwargs...)

        flow_est = outputs[1]
        source_reg = outputs[2]
    end # "timer"

    # get total runtime
    runtime = TimerOutputs.time(timer[timer.name])/10e8
    # test flow quality
    flow_names_vals_dict = asses_flow_quality(flow, flow_est)
    # test source_reg quality
    img_names_vals_dict = asses_source_reg_quality(target, source_reg)

    results = merge(flow_names_vals_dict, img_names_vals_dict)
    results["time"] = runtime

    if display
        print_timer(timer)
        println()
        print_dict(results)
    end

    if length(outputs) >= 3
        return flow_est, source_reg, timer, results, outputs[3:end]
    end

    return flow_est, source_reg, timer, results
end
