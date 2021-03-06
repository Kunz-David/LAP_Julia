image_1 = img
image_2 = imgw
filter_half_size = fhs
window_size = window_size

## function
image_size = size(image_1)
pixel_count = length(image_1)
point_count = size(points, 2)
filter_size = 2 * filter_half_size + 1

# Prepare filters:
# Calculate separable filters from basis:
sigma = (filter_half_size + 2) / 4
centered_inds = centered(-filter_half_size:filter_half_size)
gaus = gaussian(sigma, filter_size)
gaus_der = gaus .* centered_inds .* (-1)/sigma^2
gaus_der_flip = reverse(gaus_der, dims=1)

# prepare 2D filter basis for later use
basis = zeros(filter_size, filter_size, filter_num)

# dif = (points of forward image - points of backward image) for each filter
dif = zeros(image_size..., filter_num)

if filter_num == 3
    # temporary place to store filtered images
    tmp_filtered_1 = zeros(image_size)
    tmp_filtered_2 = zeros(image_size)
    tmp_filtered_1 = fill(NaN, image_size)
    tmp_filtered_2 = fill(NaN, image_size)

    # basis 1 - (gaus, gaus) both forward and backward
    kernf_1 = kernelfactors((gaus, gaus))
    basis[:, :, 1] = broadcast(*, kernf_1...)
    filt_onebyone!(tmp_filtered_1, image_1, kernf_1, filter_half_size, points)
    filt_onebyone!(tmp_filtered_2, image_2, kernf_1, filter_half_size, points)
    dif[:, :, 1] = tmp_filtered_2 - tmp_filtered_1

    # basis 2 - (gaus, gaus_der_flip) as forward, (gaus, gaus_der) as backward
    kernf_2f = kernelfactors((gaus, gaus_der_flip))
    basis[:, :, 2] = broadcast(*, kernf_2f...)
    filt_onebyone!(tmp_filtered_1, image_1, kernf_2f, filter_half_size, points)
    kernf_2b = kernelfactors((gaus, gaus_der))
    filt_onebyone!(tmp_filtered_2, image_2, kernf_2b, filter_half_size, points)
    dif[:, :, 2] = tmp_filtered_2 - tmp_filtered_1

    # basis 3 - (gaus_der_flip, gaus) as forward, (gaus_der, gaus) as backward
    kernf_3f = kernelfactors((gaus_der_flip, gaus))
    basis[:, :, 3] = broadcast(*, kernf_3f...)
    filt_onebyone!(tmp_filtered_1, image_1, kernf_3f, filter_half_size, points)
    kernf_3b = kernelfactors((gaus_der, gaus))
    filt_onebyone!(tmp_filtered_2, image_2, kernf_3b, filter_half_size, points)
    dif[:, :, 3] = tmp_filtered_2 - tmp_filtered_1
end

dif = reshape(dif, (:, filter_num))

# prepare matrices for linear system of equations
A = zeros(filter_num-1, filter_num-1, pixel_count)
b = zeros(filter_num-1, pixel_count)
A = fill(NaN, filter_num-1, filter_num-1, pixel_count)
b = fill(NaN, filter_num-1, pixel_count)

for k in 1:filter_num-1
    for l in k:filter_num-1
        @views window_sum_around_points!(A[k, l, :], dif[:, k+1] .* dif[:, l+1], image_size, window_size, points)
        A[l, k, :] = A[k, l, :]
    end
    @views window_sum_around_points!(b[k, :], dif[:, k+1] .* dif[:, 1] .* (-1), image_size, window_size, points)
end

println(length(filter(!isnan, real.(A[1, 1, :]))))
imgshow(reshape(A[1, 1, :], image_size...))
println(filter(!isnan, real.(A[1, 1, :])))
println(filter(!isnan, real.(A[2, 1, :])))
println(filter(!isnan, real.(A[2, 2, :])))
println(filter(!isnan, real.(A[1, 2, :])))

# Perform Gauss elimination on all pixels in parallel:
# coeffs will be of shape: pixel_count, filter_num-1
@views coeffs = multi_mat_div_at_points(A, b, points, image_size)
# adding ones so that all base filters have their coefficients even the first one
all_coeffs = [ones(pixel_count) coeffs]

# make a border mask. Is 0 if its in a border of filter_half_size size.
window_half_size = Int64.((window_size .- 1) ./ 2)
border_mask = parent(padarray(ones((image_size .- (2 .* window_half_size))...),
                Fill(NaN, window_half_size, window_half_size)))
@views all_coeffs = all_coeffs .* reshape(border_mask, (:, 1))

k = (-filter_half_size:filter_half_size)

# Get the displacement vector field from the filters
u1_top = zeros(image_size);
u1_bot = zeros(image_size);
u2_top = zeros(image_size);
u2_bot = zeros(image_size);

for n in 1:filter_num
    @views u1_top[:] = u1_top[:] .- sum(transpose(basis[:, :, n]) * k) .* all_coeffs[:, n];
    @views u1_bot[:] = u1_bot[:] .+ sum(basis[:, :, n]) .* all_coeffs[:, n];

    @views u2_top[:] = u2_top[:] .- sum(basis[:, :, n] * k) .* all_coeffs[:, n];
    @views u2_bot[:] = u2_bot[:] .+ sum(basis[:, :, n]) .* all_coeffs[:, n];
end

u_est = 2 .* ((im .* u1_top ./ u1_bot) .+ (u2_top ./ u2_bot));

# dont use estimations whose displacement is larger than the filter_half_size
displacement = real(u_est).^2 + imag(u_est).^2
displacement_mask = displacement .> filter_half_size^2
u_est[displacement_mask] .= NaN .+ NaN .* 1im;

any(.!isnan.(real(u_est)))

showflow(u_est)
