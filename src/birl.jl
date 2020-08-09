using DataFrames

function parse_numbers(s)
    pieces = split(s, ' ', keepempty=false)
    map(pieces) do piece
        parse(Float64, piece)
    end
end

"""
    load_H(path)

Load a homography matrix from path `path`
"""
function load_H(path)
    lines = readlines(path)
    H = Array{Float64}(undef, (3,3))
    for (k, line) in zip(1:3, lines)
        H[k,:] = parse_numbers(line)
    end
    return H
end

"""
    make_flow_from_H(H, size)

Make a flow of size `size` from the homography matrix `H`.
"""
function make_flow_from_H(H, size)
    warp_perspective(x, y) = (H[1,1]*x + H[1,2]*y + H[1,3])/(H[3,1]*x + H[3,2]*y + H[3,3]),
                             (H[2,1]*x + H[2,2]*y + H[2,3])/(H[3,1]*x + H[3,2]*y + H[3,3])
    flow = Array{Complex{Float64},2}(undef, size)
    X = ones(size[1]) * collect(range(0,1,length=size[2]))'
    Y = collect(range(0,1,length=size[1])) * ones(size[2])'

    X = ones(size[1]) * collect(1:size[2])'
    Y = collect(1:size[1]) * ones(size[2])'
    for ind in CartesianIndices(flow)
        calc = [X[ind], Y[ind]] .- warp_perspective(X[ind], Y[ind])
        flow[ind] = calc[1] + im*calc[2]
    end
    return flow
end


"""
    gen_rand_points(flow, point_count, mode; max_skew_const=4)

Generate `point_count` random indices over `flow` depending on selection `mode`.
"""
function gen_rand_points(flow, point_count, mode; max_skew_const=4)
    if mode == "Random_points"
        ret_inds = rand(CartesianIndices(flow), point_count)
    elseif mode == "Gridded"
        points_in_dim = round(Int64, sqrt(point_count))
        point_spacing = floor.(size(flow) ./ points_in_dim)

        ret_inds = Array{CartesianIndex{2}, 1}(undef, 0)
        for ind in CartesianIndices(flow)
            if all((Tuple(ind)) .% point_spacing .== 0)
                add_ind = CartesianIndex(Tuple(ind) .- floor.(Int64, point_spacing ./ 2))
                push!(ret_inds, add_ind)
            end
        end
    elseif mode == "Semigridded"
        points_in_dim = round(Int64, sqrt(point_count))
        point_spacing = floor.(size(flow) ./ points_in_dim)

        ret_inds = Array{CartesianIndex{2}, 1}(undef, 0)
        for ind in CartesianIndices(flow)
            if all((Tuple(ind)) .% point_spacing .== 0)
                skew = (rand((-floor(Int64, point_spacing[1]/max_skew_const):floor(Int64, point_spacing[1]/max_skew_const))),
                        rand((-floor(Int64, point_spacing[1]/max_skew_const):floor(Int64, point_spacing[2]/max_skew_const))))
                add_ind = CartesianIndex(Tuple(ind) .- floor.(Int64, point_spacing ./ 2) .+ skew)
                @assert add_ind in CartesianIndices(flow)

                push!(ret_inds, add_ind)
            end
        end
    end
    return ret_inds
end


function save_shift_landmarks(source_landmarks_path, flow, args_dict)
    source_landmarks_df = CSV.read(source_landmarks_path) |> DataFrame
    locations = [source_landmarks_df[:, "Y"] source_landmarks_df[:, "X"]]
    shifted_landmarks = LAP_julia.move_landmarks(locations, flow)
    shifted_landmarks_df = DataFrame(Column1 = 0:(size(locations,1)-1), X = shifted_landmarks[:,2], Y = shifted_landmarks[:,1])
    save_path = joinpath(args_dict["output_path"], args_dict["land_warped_fname"])
    CSV.write(save_path, shifted_landmarks_df)
end

function save_shift_landmarks(source_landmarks_path, flow, save_path::String)
    source_landmarks_df = CSV.read(source_landmarks_path) |> DataFrame
    locations = [source_landmarks_df[:, "Y"] source_landmarks_df[:, "X"]]
    shifted_landmarks = LAP_julia.move_landmarks(locations, flow)
    shifted_landmarks_df = DataFrame(Column1 = 0:(size(locations,1)-1), X = shifted_landmarks[:,2], Y = shifted_landmarks[:,1])
    save_path = joinpath(save_path)
    CSV.write(save_path, shifted_landmarks_df)
end

function save_landmarks(inds, save_path::String)
    points = transpose(inds_to_points(inds))
    shifted_landmarks_df = DataFrame(Column1 = 0:(size(points,1)-1), X = points[:,2], Y = points[:,1])
    save_path = joinpath(save_path)
    CSV.write(save_path, shifted_landmarks_df)
end
