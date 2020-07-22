using LAP_julia: average_dicts, print_dict

# prepare new dataframe
calculate_flow_df = DataFrame(
    index = Int[],
    mode = Symbol[],
    img_size = Int[],
    whs = Int[],
    timer = TimerOutput[],
    results = Dict{String,Float64}[]
    )

window_half_sizes = cat(collect(1:51), collect(61:10:101), collect(141:40:381), dims=1);
img_sizes = [50, 100, 200, 400, 800, 1600];
modes = [:known_value; :calculation];

##
# fill dataframe
df = calculate_flow_df
let index = 0
    for img_size in img_sizes
        img, imgw, flow = gen_init(:chess, chess_args=[25, img_size/25])

        # limit the max window half size to 1/4 the image
        whs_limit = img_size/4
        whs_modified = filter(x -> x <= whs_limit, window_half_sizes)

        for whs in whs_modified
            for mode in modes

                window = [whs * 2 + 1, whs * 2 + 1]
                timer = TimerOutput("reg alg: lap")
                # run once:
                flow_est, source_reg, timer, results = test_registration_alg(lap, img, imgw, flow, [whs, window],
                    Dict(:timer => timer, :new_feature => mode == :calculation ? false : true), timer=timer, display=true)

                index = index + 1
                println("at index: ", index, " img_size: ", img_size, " whs: ", whs)
                push!(df, Dict(:index => index,
                               :mode => mode,
                               :img_size => img_size,
                               :whs => whs,
                               :timer => timer,
                               :results => results))
            end
        end
    end
end

new_avg_dicts = []
old_avg_dicts = []

for img_size in img_sizes

    n = count(k -> k < (img_size/4), window_half_sizes)

    new_dicts = df[(df.mode .== :known_value) .& (df.img_size .== img_size), :][:, :results][1:n]
    old_dicts = df[(df.mode .== :calculation) .& (df.img_size .== img_size), :][:, :results][1:n]

    new_avg = average_dicts(new_dicts)
    old_avg = average_dicts(old_dicts)

    append!(new_avg_dicts, new_avg)
    append!(old_avg_dicts, old_avg)

    print_dict(new_avg, "new_avg, with img_size: " * string(img_size))
    print_dict(old_avg, "old_avg, with img_size: " * string(img_size))
end




img_size = 50
n = count(k -> k < (img_size/4), window_half_sizes)
new_dicts = df[(df.mode .== :known_value) .& (df.img_size .== img_size), :][:, :results][1:n]
old_dicts = df[(df.mode .== :calculation) .& (df.img_size .== img_size), :][:, :results][1:n]
