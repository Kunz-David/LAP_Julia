function [uest,coeffs,err]=optiflowFilter2(I1,I2,basis, K1, W)
% Function implements the standard Local All-Pass Filter (LAP) algorithm
% (i.e. no multi-scale). Given input images I1 and I2, and basis (the set
% of filters), the function calculates the optical flow between the two
% images. Currently designed to use grayscale images.
%
% Inputs:
%           I1      -> Input image 1, gray scale, s1 by s2 array
%           I2      -> Input image 2, gray scale, s1 by s2 array
%           basis   -> Elements of the filter basis. M by N array, where N
%                   is the number of basis elements and M = R^2 where is
%                   the filter is R by R in size.
%           K1     -> Size of filter basis
%           W      -> Size of the local window used in the algorithm.
%
% Outputs:
%           uest    -> Estimate of the optical flow as a complex number 
%           coeffs  -> Estimated coefficients for each filter basis
%           err     -> Error for each pixel: ||p(-x)*I2(x) - p(x)*I1(x)||^2
%
% Version 1.
% DATE     : 21st October 2016
% AUTHOR   : Christopher Gilliam (dr.christopher.gilliam@ieee.org) and
%            Thierry Blu
%
% Version 2.
%           - Window size can be defined independently of the filter size.
%           - Separable 1D filter implementation for large filter sizes
%               assuming N = 3.
%           - Set unreliable edge pixels to nans
% DATE     : 4th January  2017
%
% Code relates to following papers:
% 1) C. Gilliam and T. Blu, "Local All-Pass Geometric Deformations", IEEE
%   Trans. Image Processing, 2016, under review.
% 2) C. Gilliam and T. Blu, "Local All-Pass Filters for Optical Flow
%   Estimation", in Proc. IEEE Int. Conf. on Acoustics, Speech & Signal
%   Processing (ICASSP 2015), Brisbane, Australia, April 2015, pp. 1533-1537.
% 3) T. Blu, P. Moulin and C. Gilliam, "Approximation Order of the LAP
%   Optical Flow Algorithm?" in Proc. IEEE Int. Conf. Image Processing
%   (ICIP 2015), Quebec City, Canada, September 2015, pp. 48-52.
%

% Obtain image size:
[s1,s2]=size(I1);
[M,N]=size(basis);

if nargin < 4,
    % set default values of K1 and W
    K1 = (ones(2,1)*sqrt(M)-1)/2;
    W = 2.*ceil(K1)+1;
elseif nargin < 5,
    % If not specified, set window size equal to filter size:
    W = 2.*ceil(K1)+1;
end

K0 = 2*ceil(K1(1))+1;
L0 = 2*ceil(K1(2))+1;
basis=reshape(basis,[K0,L0,N]);
K=(K0-1)/2;
L=(L0-1)/2;
if round(K)~=K|round(L)~=L
    error('Basis filters are not centered.')
end

% Test size of basis filters to decided whether to use 2D  imfilter or
% break the filter into separable 1D filters (allows for faster implementation):
if K >= 64 || L  >= 64 && N == 3,
    FT_trig = 1;

    % Calculate separable filters from basis:
    sigma=(K1 + 2)/4;
    g=@(x,s)exp(-x.*x/2/s^2);
    K0=ceil([K,L]);
    Gx= g(-K0(1):K0(1),sigma(1));
    Gdx= (-K0(1):K0(1)).*Gx;
    Gx=Gx/norm(Gx,1);
    Gdx=Gdx/norm((-K0(1):K0(1)).*Gdx,1);
    Gix=Gx(end:-1:1); % inverse, although they are same
    Gdix=Gdx(end:-1:1);

    Gy= g(-K0(2):K0(2),sigma(2));
    Gdy= (-K0(2):K0(2)).*Gy;
    Gy=Gy/norm(Gy,1);
    Gdy=Gdy/norm((-K0(2):K0(2)).*Gdy,1);
    Giy=Gy(end:-1:1); % inverse, although they are same
    Gdiy=Gdy(end:-1:1);
else
    FT_trig = 0;
end

K1=W(1);         % block size first dimension
K2=W(2);              % block size second dimension

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% filtering with the basis %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
II=[];

for n=1:N
    if FT_trig == 0,
        % standard 2D filtering
        B1=basis(:,:,n);
        B2=B1(end:-1:1,end:-1:1);
        J=imfilter(I2,B2,'symmetric')-imfilter(I1,B1,'symmetric');
    else
        switch n
            case 1
                % Instead of using a 2D gaussian kernel,
                % uses the fact that a gaussian filter can be separated in 1D gaussian kernels.
                A1 = imfilter(I1, reshape(Gx,[length(Gx) 1]),'symmetric');
                A1 = imfilter(shiftdim(A1,1), reshape(Gy,[length(Gy) 1]),'symmetric').';
                B1 = imfilter(I2, reshape(Gix,[length(Gix) 1]),'symmetric');
                B1 = imfilter(shiftdim(B1,1), reshape(Giy,[length(Giy) 1]),'symmetric').';
            case 2
                A1 = imfilter(I1, reshape(Gx,[length(Gx) 1]),'symmetric');
                A1 = imfilter(shiftdim(A1,1), reshape(Gdy,[length(Gdy) 1]),'symmetric').';
                B1 = imfilter(I2, reshape(Gix,[length(Gix) 1]),'symmetric');
                B1 = imfilter(shiftdim(B1,1), reshape(Gdiy,[length(Gdiy) 1]),'symmetric').';
            case 3
                A1 = imfilter(I1, reshape(Gdx,[length(Gdx) 1 1]),'symmetric');
                A1 = imfilter(shiftdim(A1,1), reshape(Gy,[length(Gy) 1 1]),'symmetric').';
                B1 = imfilter(I2, reshape(Gdix,[length(Gdix) 1]),'symmetric');
                B1 = imfilter(shiftdim(B1,1), reshape(Giy,[length(Giy) 1]),'symmetric').';
        end
        J=B1-A1;
    end
    II=[II J(:)];
end
J=II(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% matrices needed %%%%%%%%%
%%%%%% in the linear system %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=zeros(N-1,N-1,s1*s2);
b=zeros(N-1,s1*s2);
for k=1:N-1
    for l=k:N-1
        A(k,l,:)=average(II(:,k+1).*II(:,l+1));
        A(l,k,:)=A(k,l,:);
    end
    b(k,:)=-average(II(:,k+1).*J);
end

% Initialise filter coefficients:
coeffs=zeros(s1*s2,N-1);

% Perform Gauss elimination on all pixels in parallel:
for k=1:(N-1)
    for l=(k+1):(N-1)
        c=A(l,k,:)./A(k,k,:);
        for m=(k+1):(N-1)
            A(l,m,:)=A(l,m,:)-c.*A(k,m,:);
        end
        A(l,k,:)=0;
        b(l,:)=b(l,:)-shiftdim(c,1).*b(k,:);
    end
end
for k=(N-1):-1:1
    coeffs(:,k)=shiftdim(b(k,:));
    for m=(k+1):(N-1)
        coeffs(:,k)=coeffs(:,k)-shiftdim(A(k,m,:)).*coeffs(:,m);
    end
    coeffs(:,k)=coeffs(:,k)./shiftdim(A(k,k,:));
end

coeffs=[ones(s1*s2,1) coeffs];

% Exclusion of the boundaries
p=bndindex(I1,[(K1-1)/2,(K2-1)/2]);
coeffs(p,:)=NaN;

if nargout>2
    err=zeros(s1,s2);
    err(:)=sum(coeffs.*II,2);
end

k0=(-K:K)'*ones(1,2*L+1);
l0=ones(2*K+1,1)*(-L:L);
basis=reshape(basis,[(2*K+1)*(2*L+1) N]);

% Determine displacement vector from the local all-pass filters
u1=zeros(s1,s2);
u11=zeros(s1,s2);
u2=zeros(s1,s2);
u22=zeros(s1,s2);
for n=1:N
    u1(:)=u1(:)-(basis(:,n)'*k0(:))*coeffs(:,n);
    u11(:)=u11(:)+sum(basis(:,n))*coeffs(:,n);
    u2(:)=u2(:)-(basis(:,n)'*l0(:))*coeffs(:,n);
    u22(:)=u22(:)+sum(basis(:,n))*coeffs(:,n);
end
uest=2*(u1./(u11)+i*u2./(u22));
test = real(uest).^2./(L^2) + imag(uest).^2./(K^2);
uest(find(test>1))=NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Embedded functions %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J=average(I)
K1=evalin('caller','K1');
K2=evalin('caller','K2');
s1=evalin('caller','s1');
s2=evalin('caller','s2');
FT_trig = evalin('caller','FT_trig');

if FT_trig == 1,
    F1 = ones(1,K1);
    F1 = F1./norm(F1(:),1);
    F2 = ones(1,K2);
    F2 = F2./norm(F2(:),1);

    J = imfilter(reshape(I,s1,s2), reshape(F1,[length(F1) 1]),'symmetric');
    J = imfilter(shiftdim(J,1), reshape(F2,[length(F2) 1]),'symmetric').';

else
    J=imfilter(reshape(I,s1,s2),ones(K1,K2),'symmetric')/(K1*K2);
end
J=J(:);
return

function p=bndindex(I,K)
if length(K)==1
    K=[K,K];
end
[M,N]=size(I);
m=(1:M)'*ones(1,N);
n=ones(M,1)*(1:N);
p=find(~(m>=K(1)+1&m<=M-K(1)&n>=K(2)+1&n<=N-K(2)));
return
