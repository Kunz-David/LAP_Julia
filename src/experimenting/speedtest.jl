using TimerOutputs


"""
    test_registration_alg(args; kwargs)

Test a registration function `method` by timing it and comparing its outputs to the ground truth `flow` and the `target` image.

# Arguments:
- `method`: registration function.
- `target::Image`: target image.
- `source::Image`: source image to be warped closer to `target`.
- `flow::Flow`: truth displacement flow warping `target` to `source`.

# Keyword Arguments:
- `display::Bool=true`: print debugging info and test results.
- `timer::TimerOutput=TimerOutput("blank")`: timer to for certain blocks of code in the `method`.
- `method_args=[]`: arguments for the registration function.
- `method_kwargs::Dict=Dict()`: keyword arguments for the registration function.


# Outputs:
- `flow_est`: estimated transformation flow.
- `source_reg`: `source` image registered.
- `timer`: timings of the `method` insides.
- `results`: Dict with results of the quality of the registration.
- [`(outputs)`]: other outputs of the `method` besides `flow_est` and `source_reg`.

# Example
```@example
# generate a lena image, a warped lena image and a random quadratic flow
img, imgw, flow = gen_init();
# setup timer
timer=TimerOutput("sparse pf lap");
# choose method params
method_kwargs = Dict(:timer => timer, :display => false, :max_repeats => 1, :point_count => 500, :spacing => 10)
# run and test `sparse_pflap`
flow_est, source_reg, timer, results = test_registration_alg(sparse_pflap, img, imgw, flow, [], method_kwargs, timer=timer)
```

See also: [`lap`](@ref), [`sparse_lap`](@ref),[`pflap`](@ref), [`sparse_pflap`](@ref), [`Flow`](@ref), [`time_reg_alg`](@ref)
"""
function test_registration_alg(method,
                               target::Image,
                               source::Image,
                               flow::Flow;
                               method_args=[],
                               display::Bool=true,
                               timer::TimerOutput=TimerOutput("blank"),
                               only_flow_compare::Bool=true,
                               method_kwargs::Dict=Dict(:timer => timer, :display => display))
    TimerOutputs.enable_debug_timings(LAP_julia)

    @timeit timer timer.name begin
        outputs = method(target, source, method_args...; timer=timer, method_kwargs...)

        flow_est = outputs[1]
        source_reg = outputs[2]
    end # "timer"

    # get total runtime
    runtime = TimerOutputs.time(timer[timer.name])/10e8
    # test flow quality
    flow_names_vals_dict = assess_flow_quality(flow.*(-1), flow_est)
    # test source_reg quality
    if !only_flow_compare
        img_names_vals_dict = assess_source_reg_quality(target, source_reg)
        results = merge(flow_names_vals_dict, img_names_vals_dict)
    else
        results = flow_names_vals_dict
    end

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


"""
    time_reg_alg(args; kwargs)

Run a registration method and get the time of execution and method outputs.

# Arguments:
- `method`: registration function.
- `target::Image`: target image.
- `source::Image`: source image to be warped closer to `target`.

# Keyword Arguments:
- `method_args=[]`: arguments for the registration function.
- `display::Bool=true`: print debugging info and test results.
- `timer::TimerOutput=TimerOutput("blank")`: timer to for certain blocks of code in the `method`.
- `method_kwargs::Dict=Dict(:timer => timer, :display => display)`: keyword arguments for the registration function.

# Outputs:
- `flow_est`: estimated transformation flow.
- `source_reg`: `source` image registered.
- `timer`: timings of the `method` insides.
- `time_in_secs`: execution time of the method in seconds.
- [`(outputs)`]: other outputs of the `method` besides `flow_est` and `source_reg`.


See also: [`test_registration_alg`](@ref)
"""
function time_reg_alg(method,
                      target::Image,
                      source::Image;
                      method_args=[],
                      display::Bool=true,
                      timer::TimerOutput=TimerOutput("blank"),
                      method_kwargs::Dict=Dict(:timer => timer, :display => display))

    TimerOutputs.enable_debug_timings(LAP_julia)
    # TimerOutputs.disable_debug_timings(LAP_julia)

    @timeit timer timer.name begin
        outputs = method(target, source, method_args...; timer=timer, method_kwargs...)

        flow_est = outputs[1]
        source_reg = outputs[2]
    end # "timer"

    # get total runtime
    time_in_secs = TimerOutputs.time(timer[timer.name])/10e8

    if display
        print_timer(timer)
    end

    if length(outputs) >= 3
        return flow_est, source_reg, timer, time_in_secs, outputs[3:end]
    end

    return flow_est, source_reg, timer, time_in_secs
end
