using LAP_julia
using TimerOutputs

function asses_source_reg_quality(target, source_reg; title="", display::Bool=true)
    functions = [ncc, mae, rmse, mse]
    names_short = ["ncc", "mae", "rmse", "mse"]
    vals = map(x -> x(target, source_reg), functions)
    names_vals_dict = Dict(names_short[i] => vals[i] for i in 1:4)
    return names_vals_dict
end

function asses_flow_quality(flow, flow_est; title="", display::Bool=true)
    functions = [angle_mae, angle_rmse, mae, mse]
    names_short = ["angle-mae", "angle-rmse", "mae", "rmse"]
    # names = ["angle mean absolute error", "angle root mean squared error", "mean absolute error", "root mean squared error"]
    vals = map(x -> x(flow, flow_est), functions)
    names_vals_dict = Dict(names_short[i] => vals[i] for i in 1:4)
    return names_vals_dict
end

function print_results(names_vals_dict, title="")
    max_len = maximum(length.(keys(names_vals_dict)))
    line_len = max_len + 5
    if title != ""
        println(repeat("-", line_len+10))
        println(repeat(" ", 2) * title)
        println(repeat("-", line_len+10))
    end
    lines = [repeat(" ", 2) * x * repeat(" ", line_len-length(x)) * "| " * string(round(y, digits=3)) for (x, y) in names_vals_dict]
    map(x -> println(x), lines)
end


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
    runtime = TimerOutputs.tottime(timer)/10e8
    # test flow quality
    flow_names_vals_dict = asses_flow_quality(flow, flow_est)
    # test source_reg quality
    img_names_vals_dict = asses_source_reg_quality(target, source_reg)

    results = merge(flow_names_vals_dict, img_names_vals_dict)
    results["time"] = runtime

    if display
        print_timer(timer)
        println()
        print_results(results)
    end

    return flow_est, source_reg, timer, results
end


# for reg_fun in [sparse_lap, lap]
#     s = Symbol(reg_fun)
#     println(s)
# end
