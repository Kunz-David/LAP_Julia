
function print_dict(results_dict, title="")
    max_len = maximum(length.(keys(results_dict)))
    line_len = max_len + 5
    if title != ""
        println(repeat("-", line_len+10))
        println(repeat(" ", 2) * title)
        println(repeat("-", line_len+10))
    end
    lines = [repeat(" ", 2) * x * repeat(" ", line_len-length(x)) * "| " * string(round(y, digits=3)) for (x, y) in results_dict]
    map(x -> println(x), lines)
end
