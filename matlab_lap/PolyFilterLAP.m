function [u_est,source_reg] = PolyFilterLAP(target, source)
% The function implements a Poly-Filter framework for the LAP registration 
% algorithm. Instead of downsampling the images, the framework changes the 
% size of the all-pass filters used in the LAP algorithm. The filter basis
% used in the LAP algorithm spans the derivatives of a Gaussian filter.
% Note that this implementation is for greyscale images only.
%
% Key Elements: i) Uses filter pyramid implementation
%               ii) Allow inter-scale iterations at the same filter size
%               iii) Gaussian filtering at all levels
%               iv) Able to deal with image noise
%               v) Median filtering at fine levels (amp <= 2)
%               vi) Shifted-linear interpolation at coarse levels 
%               vii) Cubic OMOMS interpolation at fine levels (amp <= 2)
%               viii) Grayscale Images
%
% USAGE : [u_est,source_reg] = PolyFilterLAP(target, source);
% INPUT : target      -> target/fixed image, gray scale, M by N array
%         source      -> source/moving image, gray scale, M by N array
%
% OUTPUT: u_est       -> 2D estimate of the displacement, M by N array
%                        (Note that a complex notation is used u = ux + j*uy)
%         source_reg  -> registered image (gray scale, M by N array), 
%                        obtained by interpolation_kernel the source image using u_est
%
% AUTHOR: Christopher Gilliam (dr.christopher.gilliam@ieee.org)
%         Thierry Blu (thierry.blu@m4x.org)
% DATE  : 4th January  2017 
%
% REFS  :
% 1) C. Gilliam and T. Blu, "Local All-Pass Geometric Deformations", 
% IEEE Transactions on Image Processing, Vol. 27 (2), pp. 1010-1025, February 2018.
% 2) C. Gilliam and T. Blu, "Local All-Pass Filters for Optical Flow 
%   Estimation", in Proc. IEEE Int. Conf. on Acoustics, Speech & Signal 
%   Processing (ICASSP 2015), Brisbane, Australia, April 2015, pp. 1533-1537.
% 3) T. Blu, P. Moulin and C. Gilliam, "Approximation Order of the LAP 
%   Optical Flow Algorithm" in Proc. IEEE Int. Conf. Image Processing 
%   (ICIP 2015), Qu?bec, Canada, September 2015, pp. 48-52.

% Filter_num -> Number of filter basis (either 3 or 6).
Filter_num = 3;

% PreFilt -> decide whether to highpass filter the images (0 = No, 1 = Yes).
PreFilt = 0;

% MedFilt -> decide whether to median filter the flow at fine scales (0 = No, 1 = Yes).
MedFilt = 0;

% NoiseTrig -> indicate whether the images are corrupted by noise (0 = No, 1 = Yes).
NoiseTrig = 0;

% Display -> Optional parameter to Display the progress of the function (0 = No, 1 = Yes).
Display = 1;

% Level_Num -> Number of layers in the filter pyramid (binary structure 
% e.g. [5,5] => {[32,32],[16,16],[8,8],[4,4],[2,2],[1,1]}). Can be either 
% a single number (thus filters are symmetric in size) or a two numbers 
% allowing asymmetric filters.
Level_Num = floor(log2(min(size(target, 1),size(target, 2))/8)+1);

if length(Level_Num) == 1,
    Level_Num(2) = Level_Num(1);
end

% Define type of interpolation used when the filter size is small (i.e.
% accurate interpolation_kernel required at fine scales)
% interpolation_kernel = 'cubicspline';
interpolation_kernel = 'cubicOMOMS';	% Cubic OMOMS interpolation, more accurate
                                        % than cubic B-spline, and same complexity

Max_I = max([target(:); source(:)]);
Min_I = min([target(:); source(:)]);

% rescale image intensities between [0, 255]:
target = (target - Min_I)./(Max_I - Min_I).*255;
source = (source - Min_I)./(Max_I - Min_I).*255;

[M,N] = size(target);
[y, x] = meshgrid(1:N, 1:M);

if NoiseTrig
    % estimate noise level in the images (assumes uniform white additive
    % Gaussian)
    nlevel1 = (estimation_noise_variance(target));       % estimate noise variance
    nlevel2 = (estimation_noise_variance(source));

    % converts noise variance to PSNR:
    PSNR_est = 10*log10(255^2/((nlevel1+nlevel2)./2));

    % Place limit on the minimum window size:
    X = max(ceil(38 - 0.5.*PSNR_est),0);
else
    X = 1;
end

[M1,N1] = size(source);
m1 = ones(M1,N1);
if (M1 < M),
    source = [source; zeros(M-M1,N1)];
    m1 = [m1; zeros(M-M1,N1)];
elseif (M1>M),
    target = [target; zeros(M1-M,N)];
end
if (N1 < N),
    source = [source, zeros(M,N-N1)];
    m1 = [m1, zeros(M,N - N1)];
elseif (N1>N),
    target = [target, zeros(M1,N1-N)];
end    

% Initialisation:
u_est = zeros(size(target));
    
% Initialise local counter:
num_level = 0;

if (2^(Level_Num(1)+1)+1) > M,
    warning(['Level Number of ', int2str(Level_Num(1)), ' results in a filter larger than the size of the input images. Level number reduced']);
    Level_Num(1) = floor(log2(M-1)-1);
end

if (2^(Level_Num(2)+1)+1) > N,
    warning(['Level Number of ', int2str(Level_Num(2)), ' results in a filter larger than the size of the input images. Level number reduced']);
    Level_Num(2) = floor(log2(N-1)-1);
end

LN = max(Level_Num)+1;
amp_array(1,:) = 2.^(linspace(Level_Num(1),0,floor(LN)));
amp_array(2,:) = 2.^(linspace(Level_Num(2),0,floor(LN)));

index = find(sum(logical(amp_array == 2)) >= 1,1);

if isempty(index),
    index = max([1, floor(LN)-1]);
end

% threshold for iterations at a given scale:
thres = 0.3.*ones(1,length(amp_array));

% Maximum number of iterations at a given scale:
max_it = 3;


% Start estimating the displacement
if Display,h=waitbar(0);end
for l = 1:LN,
    num_level = num_level + 1;
    if Display
        waitbar(num_level/length(amp_array),h,['Level ',...
            int2str(num_level), '/',int2str(length(amp_array)),...
            ' (filter size = [', int2str(2*amp_array(1,l)+1),'x',int2str(2*amp_array(2,l)+1),'])'])
    end
    
    
    trig = 1;
    count = 1;
    if l == 1,
            % initialise parameters for the algorithm
            source_reg = source;
            PSNR_pre = 0;
    else
        if PreFilt
            % Pre-filter images and calculate the PSNR between images for
            % the new filter size before alignment at that scale
            target_prefiltered = Prefilter_Procedure(target,  amp_array(:,l));
            source_prefiltered = Prefilter_Procedure(source_reg,  amp_array(:,l));
            source_prefiltered(mask1 == 0) = target_prefiltered(mask1 == 0);
            PSNR_pre = CG_PSNR(target_prefiltered,source_prefiltered);            
        else
            % Calculate the PSNR between the images for the new filter size 
            % before alignment at that scale
            PSNR_pre = CG_PSNR(target,source_reg);
        end
    end
    
    % Inter-scale iterations with the same filter size:
    while trig == 1,

        % define scale factor
        amp_size = amp_array(:,l);
        
        % Load Filter Basis:
        Basis_Set = loadbasis_asymm(amp_size);
        Basis_Set = Basis_Set(:,1:Filter_num);
        if l == 1 && count == 1,
            if PreFilt
                % Pre-filter images using Gaussian filter
                target1 = Prefilter_Procedure(target, amp_size);
                source1 = Prefilter_Procedure(source, amp_size);
            else
                target1 = target;
                source1 = source;
            end
        else
            if PreFilt
                % Pre-filter images using Gaussian filter
                target1 = Prefilter_Procedure(target, amp_size);
                source1 = Prefilter_Procedure(source_reg, amp_size);
                source1(mask1 == 0) = target1(mask1 == 0);
            else
                source1 = source_reg;
            end
        end
                
        % calculate local window size for LAP:
        wind(1) = max([2*ceil(amp_size(1))+1,2*X+1]);
        wind(2) = max([2*ceil(amp_size(2))+1,2*X+1]);
       
        % Using basis functions estimate displacement field on a local scale:
        [uest_Orig, ~] = optiflowFilter2(target1, source1, Basis_Set, amp_size, wind);

        % Clean displacement field:
        % Stage 1: Remove edge pixels from the displacement field and replace with simple 
        % padding.
        index1 = (round(wind(1)./2)+1):(M-round(wind(1)./2));
        index2 = (round(wind(2)./2)+1):(N-round(wind(2)./2));
        uest_Clean = padarray(uest_Orig(index1,index2),[index1(1)-1,index2(1)-1],'replicate','both');
        if sum(isnan(uest_Clean(:))) < numel(uest_Clean),
            % Stage 2: Remove nan's in the displacement field using inpainting:
            [uest_Clean, ~] = cleanOF2(uest_Clean, not(isnan(uest_Clean)));
        else
            uest_Clean = zeros(M,N);
        end
            
        % Smooth flow estimation:
        uest_Clean = Cleaning_Procedure_Mean(uest_Clean, max([amp_size,X*ones(2,1)],[],2));

        % Add estimated flow to the previous flow estimate: 
        u_holder = u_est + uest_Clean;
                        
        if l < index,
            % Fast interpolation_kernel using shiftedlinear interpolation (accuracy
            % equivalent to Keys interpolation)
            source_reg = ShiftedLinear_Interp_2D(-real(u_holder), -imag(u_holder), source);
        else
            % Accurate interpolation using cubic OMOMS
            source_reg = imshift(source,-(u_holder),interpolation_kernel);
        end
        
        % generate mask to find pixels of the target image that are seen in the source image:
        mask = (logical( (y+imag(u_holder)) <= N1 & (y+imag(u_holder)) >= 1).*logical( (x+real(u_holder)) <= M1 & (x+real(u_holder)) >= 1));

        if PreFilt
            % Pre-filter images:
            target_prefiltered = Prefilter_Procedure(target, amp_size);
            source_prefiltered = Prefilter_Procedure(source_reg, amp_size);
            
            % Pixels that are not obtained by interpolation of the source
            % image are set to the target image values
            source_prefiltered(~mask) = target_prefiltered(~mask);
        else
            % Pixels that are not obtained by interpolation of the source
            % image are set to the target image values
            source_reg(~mask) = target(~mask);
            target_prefiltered = target;
            source_prefiltered = source_reg;
        end

        % Calculate PSNR for the images when the current flow has been used
        % to warp image 2 closer to image1:
        PSNR_test = CG_PSNR(target_prefiltered,source_prefiltered);
        
        if PSNR_test - PSNR_pre > thres(l),
            % if the PSNR value is above a threshold then repeat flow
            % estimation at this scale:
        
            % Update global estimate of flow:
            u_est = u_holder;
            
            % Update parameters for the algorithm
            mask1 = mask;
                       
            PSNR_pre = PSNR_test;
            count = count+1;
            if count > max_it,
                % stop inter-scale iterations at the current scale if count is greater
                % than the maximum number of iterations allowed
                trig = 0;
            end
        else
            % if the PSNR value is below a threshold then stop inter-scale iterations 
            % and reduce the filter size:
            trig = 0;
            if PSNR_test - PSNR_pre > 0,
                u_est = u_holder;
                mask1 = mask;
                PSNR_pre = PSNR_test;
            end     
        end           
   
    end
    
    if l >= index && MedFilt
        % Refinement of the displacement field at fine scales
        u_est = Cleaning_Procedure_Med(u_est);
    end
end
if Display,close(h),end
source_reg(~mask1)=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Embedded functions %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u_out = Cleaning_Procedure_Mean(u_in, amp_size)
% function cleans the estimate of the displacement field by smoothing 
% using a Gaussian filter

% Define sigma value for Gaussian filter:
R1 = round(2*amp_size(1));
R2 = round(2*amp_size(2));

% Obtain filters:
l1 = -R1:R1;
l1 = l1';
k1 = (-R2:R2);
Gauss_Filter1 = exp(-l1.*l1/2/R1^2);
Gauss_Filter1 = Gauss_Filter1./norm(Gauss_Filter1(:),1);
Gauss_Filter2 = exp(-k1.*k1/2/R2^2);
Gauss_Filter2 = Gauss_Filter2./norm(Gauss_Filter2(:),1);

% Smooth input displacement field
u_out = imfilter(u_in, reshape(Gauss_Filter1,[length(Gauss_Filter1) 1]),'symmetric');
u_out = imfilter(shiftdim(u_out,1), reshape(Gauss_Filter2,[length(Gauss_Filter2) 1]),'symmetric').';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u_out = Cleaning_Procedure_Med(u_in)
% function cleans the estimate of the displacement field, by smoothing 
% using median filters at two scales.

B1 = 15;
B2 = 5;
        
% Two part median filtering:
% 1. fine scale
u_out = medfilt2(real(u_in),[B2,B2], 'symmetric') + 1i.*medfilt2(imag(u_in),[B2,B2], 'symmetric');
% 2. coarse scale
u_out = medfilt2(real(u_out),[B1, B1], 'symmetric') + 1i.*medfilt2(imag(u_out),[B1, B1], 'symmetric');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I_out = Prefilter_Procedure(I_in, amp_size)
% function pre-filters the input image using a high-pass filter. The
% high-pass filter is 1 - Gaussian function, where the Gaussian function is
% the first basis filter. Note that the pre-filtering is performed in a
% separable manner.

% calculate the sigma values for the Gaussian function:
sigma=(amp_size + 2)/4;

% define 1D Gaussian  functions: 
g=@(x,s)exp(-x.*x/2/s^2);

K0=ceil(amp_size);
Gx= g(-K0(1):K0(1),sigma(1)); 
Gx=Gx/norm(Gx,1); 
Gy= g(-K0(2):K0(2),sigma(2)); 
Gy=Gy/norm(Gy,1); 

I_test = imfilter(I_in, reshape(Gx,[length(Gx) 1]),'symmetric');
I_out = I_in -  imfilter(shiftdim(I_test,1), reshape(Gy,[length(Gy) 1]),'symmetric').';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
