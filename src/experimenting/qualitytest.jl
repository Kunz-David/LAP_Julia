using TimerOutputs, BenchmarkTools


"""
    assess_source_reg_quality(target, source_reg; title="", display::Bool=true)

Run a few tests on comparing `target` and `source_reg`.

Tests: [`ncc`, `mae`, `rmse`, `mse`]
"""
function assess_source_reg_quality(target, source_reg; title="", display::Bool=true)
    functions = [ncc, mae, rmse, mse]
    names_short = ["ncc", "mae", "rmse", "mse"]
    vals = map(x -> x(target, source_reg), functions)
    names_vals_dict = Dict(names_short[i] => vals[i] for i in 1:4)
    return names_vals_dict
end

"""
    assess_flow_quality(flow, flow_est; title="", display::Bool=true)

Run a few tests on comparing `flow_est` and `flow`.

Tests: [`angle_mae`, `angle_rmse`, `mae`, `mse`]
"""
function assess_flow_quality(flow, flow_est; title="", display::Bool=true)
    functions = [angle_mae, angle_rmse, mae, mse]
    names_short = ["angle-mae", "angle-rmse", "flow_mae", "flow_rmse"]
    # names = ["angle mean absolute error", "angle root mean squared error", "mean absolute error", "root mean squared error"]
    vals = map(x -> x(flow, flow_est), functions)
    names_vals_dict = Dict(names_short[i] => vals[i] for i in 1:4)
    return names_vals_dict
end

# function average_dicts(dicts)
#     avg_dict = Dict()
#     bad_count = 0
#     for dict in dicts
#         for key in keys(dict)
#             if dict[key] == NaN
#                 bad_count += 1
#                 break;
#             end
#             if !haskey(avg_dict, key)
#                 avg_dict[key] = dict[key]
#             else
#                 avg_dict[key] += dict[key]
#             end
#         end
#     end
#     for key in keys(avg_dict)
#         avg_dict[key] = avg_dict[key]/(length(dicts)-bad_count)
#     end
#     return avg_dict
# end

"""
    fun_on_dict_values(dicts, fun)

Run function `fun` on all same keys of all dictionaries `dicts`.

```jldoctest
using BenchmarkTools
fun_on_dict_values((Dict(:a => 1), Dict(:a => 12)), mean)

# output

Dict{Any,Any} with 1 entry:
  :a => 6.5

```
"""
function fun_on_dict_values(dicts, fun)
    out_dict = Dict()

    for key in keys(dicts[1])
        out_dict[key] = fun([dict[key] for dict in dicts])
    end
    return out_dict
end


fun_on_dict_values((Dict(:a => 1), Dict(:a => 12)), mean)

"""
    compare_dicts(old_dict, new_dict)

Compare two dicts by dividing their corresponding values, return the comparation dict.

Note: (old / new)
"""
function compare_dicts(old_dict, new_dict)
    comparation_dict = Dict()
    for key in keys(old_dict)
        if haskey(new_dict, key)
            comparation_dict[key] = old_dict[key]/new_dict[key]
        end
    end
    return comparation_dict
end
