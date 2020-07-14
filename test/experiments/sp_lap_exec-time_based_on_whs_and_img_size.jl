
using LAP_julia, TimerOutputs, DataFrames, TableView, JLD2, FileIO
##
sp_lap_win_sum = DataFrame(
    index = Int[],
    reg_fun = Symbol[],
    img_size = Int[],
    whs = Int[],
    timer = TimerOutput[],
    results = Dict{String,Float64}[]
    #Benchmark = BenchmarkTools.Trial[],
    #flow = Matrix[],
    #flow_est = Matrix[]
    )

window_half_sizes = cat(collect(1:51), collect(61:10:101), collect(141:40:381), dims=1);
img_sizes = [50, 100, 200, 400, 800, 1600];

##
df = sp_lap_win_sum
let index = 0
    for img_size in img_sizes
        img, imgw, flow = gen_init(:chess, chess_args=[25, img_size/25])

        # limit the max window half size to 1/4 the image
        whs_limit = img_size/4
        whs_modified = map(x -> x <= whs_limit ? x : 0, window_half_sizes)

        for whs in whs_modified
            for reg_fun in [sparse_lap_win_sum1, sparse_lap]
                if whs != 0
                    window = [whs * 2 + 1, whs * 2 + 1]
                    # window sum 1
                    timer = TimerOutput("reg alg: sp lap")
                    flow_est, source_reg, timer, results = test_registration_alg(reg_fun, img, imgw, flow, [whs, window], Dict(:timer => timer), timer=timer, display=false)
                else
                    timer = TimerOutput("blank")
                    results = Dict()
                end

                index = index + 1
                println("at index: ", index)
                push!(df, Dict(:index => index,
                               :reg_fun => Symbol(reg_fun),
                               :img_size => img_size,
                               :whs => whs,
                               :timer => timer,
                               :results => results))
            end
        end
    end
end

@save "sp_lap_win_sum_df.jld2" df



df

win_sum1_df = df[(df.reg_fun .== :sparse_lap_win_sum1), :]
win_sum3_df = df[(df.reg_fun .== :sparse_lap), :]

img_size = 50

sections = ["reg alg: sp lap", "sparse lap", "prepare A and b"]
sections = ["reg alg: sp lap", "sparse lap", "filtering"]

fig, axs = subplots(6, sharex=true, figsize=(10,10))
for (img_size, ax) in zip(img_sizes, axs)
    ax.set_title("Image Size: " * string(img_size))

    n = count(k -> k < (img_size/4), window_half_sizes)
    timers1 = df[(df.reg_fun .== :sparse_lap_win_sum1) .& (df.img_size .== img_size), :][:, :timer][1:n]
    timers3 = df[(df.reg_fun .== :sparse_lap) .& (df.img_size .== img_size), :][:, :timer][1:n]

    y1 = map(x -> TimerOutputs.time(get_timer(x, sections))/10e8, timers1)
    y3 = map(x -> TimerOutputs.time(get_timer(x, sections))/10e8, timers3)
    x = filter(x -> x < (img_size/4), window_half_sizes)


    xscale("log")
    yscale("log")
    subplots_adjust(top=1.5)
    ax.plot(x, y1)
    ax.plot(x, y3)
    xlabel("whs [pixels]")
    ax.set_ylabel("time [s]")
end
gcf()


df[741, :][:timer]

showtable(df)

TimerOutputs.time(get_timer(win_sum3_df[:, :timer][1], ["reg alg: sp lap", "sparse lap", "prepare A and b"]))
TimerOutputs.tottime(get_timer(win_sum3_df[:, :timer][1], ["reg alg: sp lap", "sparse lap", "prepare A and b"]))
TimerOutputs.tottime(get_timer(win_sum3_df[:, :timer][1], ["reg alg: sp lap", "interpolate flow"]))
TimerOutputs.time(win_sum3_df[:, :timer][1]["reg alg: sp lap"]["interpolate flow"])


print_timer(win_sum3_df[:, :timer][1])


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
