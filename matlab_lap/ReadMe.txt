Folder contains the Matlab code for the multi-scale implementation of the LAP algorithm. The multi-scale framework is based on filter pyramids rather than image pyramids (i.e. no downsampling of the input images).

Key Elements:

*****************************************
1) PolyFilterLAP.m
*****************************************

Standard implementation of the multi-scale LAP with the following key points:
                i) Uses filter pyramid implementation
		ii) Allows inter-scale iterations at the same filter size
		iii) Gaussian filtering at all levels
		iv) Able to deal with image noise
		v) Median filtering at fine levels (amp <= 2)
		vi) Shifted-linear interpolation at coarse levels 
		vii) Cubic OMOMS interpolation at fine levels (amp <= 2)
		viii) Grayscale Images

Example of using code:
	% Obtain input images (gralescale):
	target = double(rgb2gray(imread(-Target/Fixed Image-)));
	source = double(rgb2gray(imread(-Source/Moving Image-)));

	% Run LAP:
	[u_est,source_reg] = PolyFilterLAP(target, source);

	Note: The flow is represented using a complex number. Using the following to convert to the Middlebury format:
	uv(:,:,1) = imag(u_est);
	uv(:,:,2) = real(u_est);

*****************************************

AUTHORS : Christopher Gilliam (dr.christopher.gilliam@ieee.org)
          Thierry Blu (thierry.blu@m4x.org)
DATE    : 4th January  2017 

*****************************************
Code relates to following papers:
*****************************************

1) C. Gilliam and T. Blu, "Local All-Pass Geometric Deformations", IEEE Transactions on Image Processing, Vol. 27 (2), pp. 1010-1025, February 2018.
2) C. Gilliam and T. Blu, "Local All-Pass Filters for Optical Flow Estimation", in Proc. IEEE Int. Conf. on Acoustics, Speech & Signal Processing (ICASSP 2015), Brisbane, Australia, April 2015, pp. 1533-1537.
3) T. Blu, P. Moulin and C. Gilliam, "Approximation Order of the LAP Optical Flow Algorithm" in Proc. IEEE Int. Conf. Image Processing (ICIP 2015), Qu?bec, Canada, September 2015, pp. 48-52.
