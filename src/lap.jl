using LinearAlgebra

# Function implements the standard Local All-Pass Filter (LAP) algorithm
# (i.e. no multi-scale). Given input images I1 and I2, and basis (the set
# of filters), the function calculates the optical flow between the two
# images. Currently designed to use grayscale images.
#
# Inputs:
#           I1      -> Input image 1, gray scale, s1 by s2 array
#           I2      -> Input image 2, gray scale, s1 by s2 array
#           basis   -> Elements of the filter basis. M by N array, where N
#                   is the number of basis elements and M = R^2 where is
#                   the filter is R by R in size.
#           K1     -> Size of filter basis
#           W      -> Size of the local window used in the algorithm.
#
# Outputs:
#           uest    -> Estimate of the optical flow as a complex number
#           coeffs  -> Estimated coefficients for each filter basis
#           err     -> Error for each pixel: ||p(-x)*I2(x) - p(x)*I1(x)||^2
#

"""
    single_lap(image1, image2, base_filters, )

# input:
- `image1` ... grayscale image 1
- `image2` ... grayscale image 2
- `n_of_filters` ... ?? IS IT EVEN NEEDED HERE? NO, WILL CALL A FUNCTION
- `window` ... size of local window
- `filter_size` ... filter

"""
function single_lap(image1, image2, n_of_filters, filter_size, window)

    image_size = size(image1)
    filter_half_size = (filter_size - 1) / 2


    # check if its better to filter with 2 1D filters
    if n_of_filters == 3 && filter_size >= 129
        # Calculate separable filters from basis:
        sigma = (filter_half_size + 2) / 4

        gaus(x,s) = exp(-x.*x/2/s^2)
        # K0=ceil([K,L]); --> this is filter half size
        Gx = gaus(-filter_half_size:1:filter_half_size, sigma);
        Gdx = (-filter_half_size:1:filter_half_size) .* Gx;
        Gx=Gx ./ sum(Gx)
        Gdx=Gdx ./ norm((-K0(1):K0(1)).*Gdx,1); #--> does this need a normalization?
        # Gdix=Gdx(end:-1:1);
        #
        # % Gx is a gausian
        # % Gdx is a derivative of a gausian
        #
        # Gy= g(-K0(2):K0(2),sigma(2));
        # Gdy= (-K0(2):K0(2)).*Gy;
        # Gy=Gy/norm(Gy,1);
        # Gdy=Gdy/norm((-K0(2):K0(2)).*Gdy,1);
        # Gdiy=Gdy(end:-1:1);

    end


end
