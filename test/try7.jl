
using ArgParse, CSV, TimerOutputs, FileIO, Colors, DataFrames
using LAP_julia

println("running file...")

function save_string_to_file(str, file_path)
    open(file_path, "w") do io
        write(io, str)
    end
end

function load_image_gray(img_path)
    img = Float64.(Gray.(load(img_path)))
end

function save_shift_landmarks(source_landmarks_path, flow, args_dict)
    source_landmarks_df = CSV.read(source_landmarks_path) |> DataFrame
    locations = [source_landmarks_df[:, "Y"] source_landmarks_df[:, "X"]]
    shifted_landmarks = LAP_julia.move_landmarks(locations, flow)
    shifted_landmarks_df = DataFrame(Column1 = 0:(size(locations,1)-1), X = shifted_landmarks[:,2], Y = shifted_landmarks[:,1])
    save_path = joinpath(args_dict["output_path"], args_dict["land_warped_fname"])
    CSV.write(save_path, shifted_landmark_df)
end

function save_time(time_in_secs, args_dict)
    save_string_to_file(string(time_in_secs), joinpath(args_dict["output_path"], args_dict["time_fname"]))
end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--opt1"
            help = "an option with an argument"
        "--time_fname"
            help = "name of the file to save execution time"
            arg_type = String
            default = "exec_time_in_secs.txt"
        "--land_warped_fname"
            help = "name of the file to save warped landmarks"
            arg_type = String
            default = "warped_landmarks.csv"
        "--diag_pix"
            help = "an option without argument, i.e. a flag"
            arg_type = Int
            default = 500
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
    println("******ARGS:********")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
        args_dict[arg] = val
    end
    args_dict["reg_alg"] = getfield(LAP_julia, Symbol(args_dict["reg_alg"]))

    target = load_image_gray(args_dict["img_path"])
    source = load_image_gray(args_dict["imgw_path"])

    resized_target, target_resize_ratio = resize_to_diag_size(target, args_dict["diag_pix"])
    resized_source, source_resize_ratio = resize_to_diag_size(source, args_dict["diag_pix"])

    ready_target, ready_source = pad_images(resized_target, resized_source)

    time_in_secs = 0
    # Take the 3 run of the reg alg, allowing for julia compilation.
    for _ in 1:3
        timer = TimerOutput(string(args_dict["reg_alg"]))
        flow_est, source_reg, timer, time_in_secs = time_reg_alg(args_dict["reg_alg"],
                                                                 ready_target,
                                                                 ready_source,
                                                                 timer=timer,
                                                                 method_kwargs=Dict(:display => false, :timer => timer))
        println()

    end

    # save time
    save_time(time_in_secs, args_dict)

    # rescale flow_est
    flow_est = imresize(flow_est.*source_resize_ratio, ratio=source_resize_ratio)

    # move landmarks
    save_shift_landmarks(args_dict["moving_landmark_path"], flow_est, args_dict)
end

# img = load_image_gray("/Users/MrTrololord/Google_Drive/cvut/bakalarka/BIRL/data-images/images/artificial_reference.jpg")
# imgshow(img)

main()
