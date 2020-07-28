A = [1, 2 , 5, 6]
using Interact

@manipulate for a in slider(A, value=3, label="a")
    println(a)
end

using JLD2



@load "./test/experiments/df_point_count_save.jld2" df

df

point_counts = [50, 200, 350, 500, 650, 800, 1100, 1400]
spacings = [5, 7, 10, 14, 18, 22, 26, 30, 36, 40, 50, 60]
diag_sizes = [400, 550, 700]

sections = ["reg alg: sp lap"]
@manipulate for k in slider(1:length(spacings), value=1, label="spacing")
    spacing = spacings[k]
    fig, axs = subplots(length(diag_sizes), sharex=true, figsize=(10,10))
    for (diag_size, ax) in zip(diag_sizes, axs)
        ax.set_title("Diagonal Size: " * string(diag_size))

        local_df = df[((df.diag_size .== diag_size) .& (df.spacing .== spacing)), :]

        # get times
        timers = local_df[:, :timer]
        avg_times = map(x -> get_avg_time(get_timer(x, sections)), timers)
        println(size(avg_times))
        println(size(point_counts))
        println("--------")

        # get flow mae
        median_dicts = local_df[:, :median_results]
        median_flow_mae = map(dict -> dict["flow_mae"], median_dicts)

        # xscale("log")
        # yscale("log")
        subplots_adjust(top=1.5)
        ax.plot(point_counts, avg_times)
        # ax.plot(point_counts, median_flow_mae)
        xlabel("point count")
        ax.set_ylabel("time [s] & mae [pixel]")
    end
    gcf()
end


@manipulate for k in slider(1:length(spacings), value=1, label="spacing")
    spacing = spacings[k]
    fig, axs = subplots(3, sharex=true, figsize=(10,10), gridspec_kw=Dict("height_ratios" => [1, 1, 2]))

    axs[1].set_title("Times at spacing: " * string(spacing))
    axs[2].set_title("MAE at spacing: " * string(spacing))
    for diag_size in diag_sizes

        local_df = df[((df.diag_size .== diag_size) .& (df.spacing .== spacing)), :]

        local_point_counts = local_df[:, :point_count]

        # get times
        timers = local_df[:, :timer]
        avg_times = map(x -> get_avg_time(get_timer(x, sections)), timers)

        subplots_adjust(top=0.9)
        axs[1].plot(local_point_counts, avg_times)
        xlabel("point count")
        axs[1].set_ylabel("time [s]")

        # get flow mae
        median_dicts = local_df[:, :median_results]
        median_flow_mae = map(dict -> dict["flow_mae"], median_dicts)

        axs[2].plot(local_point_counts, median_flow_mae)
        xlabel("point count")
        axs[2].set_ylabel("mae [pixel]")

        # get points found
        points_found = local_df[:, :points_found]
        axs[3].plot(local_point_counts, points_found)
        xlabel("point count")
        axs[3].set_ylabel("points found")
        axs[3].axis("scaled")
    end
    gcf()
end



diag_sizes = [400, 550, 700];



get_avg_time(df[((df.diag_size .== 400) .& (df.spacing .== 5)), :][:, :timer][1]["reg alg: sp lap"])






##
function get_avg_time(timer)
    return (TimerOutputs.time(timer)/TimerOutputs.ncalls(timer))/10e8
end

function get_timer(timer, sections)
    private_sections = copy(sections)
    try
        timer[private_sections[1]]
    catch e
        println(e)
        return TimerOutput()
    end
    out = timer
    while private_sections != []
        out = out[popfirst!(private_sections)]
    end
    return out
end
