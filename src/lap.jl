
using Printf

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

function single_lap(image1, image2, base_filters, )
