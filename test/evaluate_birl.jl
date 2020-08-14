
using CSV, DataFrames


base_path = "/Volumes/davidkunz/Documents/bakalarka/mri_head_homography_results_all"
method_folder_names = filter(x -> !occursin(".", x), readdir(base_path))

latexify(df_summary, env=:table, latex=false)

get_mean(name) = round(Float64(results.mean[results.Column1 .== name][1]), sigdigits=3)



function get_method_summary(base_path, method_folder_name)

    method_folder = joinpath(base_path, method_folder_name)
    results_file_path = joinpath(method_folder, "results-summary.csv")
    results = CSV.read(results_file_path) |> DataFrame
    get_mean(name) = Float64(results.mean[results.Column1 .== name][1])
    # rownames = collect(results["Column1"])

    method_name = method_folder_name[findfirst(r"[a-zA-Z]+", method_folder_name)]
    method_name = method_name[1:2] == "Bm" ? method_name[3:end] : method_name
    if method_name == "Sparse"
        open(joinpath(base_path, method_folder_name, "results-summary.txt")) do f
            for i in eachline(f)
                if occursin(r"algorithm", i)
                    method_name = reverse(reverse(i)[findfirst(r"[a-zA-Z_]+", reverse(i))])
                    # @info i method_name reverse(i)
                end
            end
        end
    end
    if method_name == "pflap"
        found = false
        open(joinpath(base_path, method_folder_name, "0", "registration.log")) do f
            for i in eachline(f)
                if occursin(":prefilter => true", i)
                    found = true
                    # @info i method_name reverse(i)
                end
            end
        end
        resized = false
        open(joinpath(base_path, method_folder_name, "results-summary.txt")) do f
            for i in eachline(f)
                if occursin(r"diag_pix", i)
                    resized = true
                    global reduced_size = string(i[findfirst(r"[0-9]+", i)])
                end
            end
        end
        if found
            if resized
                method_name = "\$\\text{PF-LAP}_{\\text{prefilt}}\$\\textsuperscript{*}"
            else
                method_name = "\$\\text{PF-LAP}_{\\text{prefilt}}\$"
            end
        else
            if resized
                method_name = "\$\\text{PF-LAP}\$\\textsuperscript{*}"
            else
                method_name = "\$\\text{PF-LAP}\$"
            end
        end
    elseif method_name == "sparse_pflap"
        resized = false
        open(joinpath(base_path, method_folder_name, "results-summary.txt")) do f
            for i in eachline(f)
                if occursin(r"diag_pix", i)
                    resized = true
                    global reduced_size = string(i[findfirst(r"[0-9]+", i)])
                end
            end
        end
        if resized
            method_name = "\$\\text{Sparse PF-LAP}\\textsuperscript{*}\$"
        else
            method_name = "\$\\text{Sparse PF-LAP}\$"
        end
    elseif method_name == "sparse_pflap_psnr"
        resized = false
        open(joinpath(base_path, method_folder_name, "results-summary.txt")) do f
            for i in eachline(f)
                if occursin(r"diag_pix", i)
                    resized = true
                    global reduced_size = string(i[findfirst(r"[0-9]+", i)])
                end
            end
        end
        if resized
            method_name = "\$\\text{Sparse PF-LAP}_{\\mathrm{PSNR}}\\textsuperscript{*}\$"
        else
            method_name = "\$\\text{Sparse PF-LAP}_{\\mathrm{PSNR}}\$"
        end
    else
        method_name = "\$\\text{$method_name}\$"
    end

    exec_time = "Execution time [minutes]"
    exec_time_secs = "Execution time [seconds]"
    tre_mean = "TRE Mean"
    tre_median = "TRE Median"
    diag_pix = "Image diagonal [pixels]"

    wanted_values = [tre_mean, tre_median, exec_time, diag_pix]

    df_summary = DataFrame(Symbol("Method") => method_name,
                           Symbol("\$\\mathrm{TRE}_{\\text{Mean}}\$") => round(get_mean(tre_mean), sigdigits=3),
                           Symbol("\$\\mathrm{TRE}_{\\text{Med}}\$") => round(get_mean(tre_median), sigdigits=3),
                           Symbol("Time") => round(get_mean(exec_time)*60, sigdigits=3)
                           # Symbol(diag_pix) => get_mean(diag_pix))
                           )
    return df_summary
end

function get_df_all_methods(base_path)
    method_folder_names = filter(x -> !occursin(".", x), readdir(base_path))

    df = get_method_summary(base_path, method_folder_names[1])
    delete!(df, 1:size(df,1))

    for folder_name in method_folder_names
        method_df = get_method_summary(base_path, folder_name)
        # @info eachcol(method_df)
        push!(df, DataFrameRow(method_df, 1))
    end
    return df
end



base_path = "/Volumes/davidkunz/Documents/bakalarka/ox_affine_separate_results/bikes"
base_path = "/Volumes/davidkunz/Documents/bakalarka/ox_affine_separate_results/ubc"
df_out = get_df_all_methods(base_path)

describe(df_out)



base_path = "/Volumes/davidkunz/Documents/bakalarka/mri_head_homography_results_all/"
df_out = get_df_all_methods(base_path)

describe(df_out)

using Latexify

latexify(df_out, env=:table, latex=false)

Dict(eachrow(get_df_summary(base_path, method_folder_names[4])))


## ox affine tables

folders = ["bikes", "leuven", "trees", "ubc"]
base_path = "/Volumes/davidkunz/Documents/bakalarka/ox_affine_separate_results"

base_path = "/Volumes/davidkunz/Documents/bakalarka/"
folders = ["mri_head_homography_results_all"]


for folder in folders
    path = joinpath(base_path, folder)
    df_out = get_df_all_methods(path)
    @info folder
    println(latexify(df_out, env=:table, latex=false))
end
