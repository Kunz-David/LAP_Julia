
using ArgParse, CSV, TimerOutputs, FileIO, Colors
using LAP_julia
println("hello world")


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--opt1"
            help = "an option with an argument"
        "--time_fname", "-o"
            help = "name of the file to save execution time"
            arg_type = String
            default = "exec_time_in_secs.txt"
        "--flag1"
            help = "an option without argument, i.e. a flag"
            action = :store_true
        "img_path"
            help = "target image path"
            required = true
        "imgw_path"
            help = "source image path"
            required = true
        "moving_landmark_path"
            help = "source landmark path"
            required = true
        "output_path"
            help = "output path"
            required = true
        "reg_alg"
            help = "registration algorithm"
            required = true
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println("Running with the arguments:")
    args_dict = Dict()
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
        args_dict[arg] = val
    end
    args_dict["reg_alg"] = getfield(LAP_julia, Symbol(args_dict["reg_alg"]))

    target = load_image_gray(args_dict["img_path"])
    source = load_image_gray(args_dict["imgw_path"])

    source_landmarks = CSV.read(args_dict["moving_landmark_path"])

    time_in_secs = 0

    # Take the 3 run of the reg alg, allowing for julia compilation.
    for _ in 1:3
        timer = TimerOutput(string(args_dict["reg_alg"]))
        flow_est, source_reg, timer, time_in_secs = time_reg_alg(args_dict["reg_alg"],
                                                                 target,
                                                                 source,
                                                                 timer=timer,
                                                                 method_kwargs=Dict(:display => false, :timer => timer))
        println()

    end

    # save time
    save_string_to_file(string(time_in_secs), joinpath(args_dict["output_path"], args_dict["time_fname"]))

    # move landmarks

end

function save_string_to_file(str, file_path)
    open(file_path, "w") do io
        write(io, str)
    end
end


function load_image_gray(img_path)
    img = Float64.(Gray.(load(img_path)))
end

# img = load_image_gray("/Users/MrTrololord/Google_Drive/cvut/bakalarka/BIRL/data-images/images/artificial_reference.jpg")
# imgshow(img)


main()
