function [u0,mask]=cleanOF2(u,mask,display)
% USAGE: [u0,mask]=cleanOF(u,mask)
% mask is an image with same size as u equal to
%       * 0 for the pixels of u to be inpainted
%       * 1 for the pixels of u to be kept as is
%       * 2 for edge pixels
if nargin<2
    mask=masking(u,'Laplacian');
    display = 0;
elseif nargin <3,
    display = 0;
end
mask=double(mask);
u(find(mask~=1))=0;

% Diffusion

% h=[0 1 0;1 0 1;0 1 0];h=conv2(h,h);
b=1;a=b/4;h=[a b a;b 0 b;a b a];
% h=ones(3,3);h(2,2)=0;
% K0=1;s=K0/2;g=@(x)exp(-x.*x/2/s^2);k=-K0:K0;k=k'*ones(1,2*K0+1);l=k';h=g(k).*g(l);h(K0+1,K0+1)=0;

isedge=(mask==2);
noedge=find(~isedge);

u0=u;
encore=(sum(mask(:))>=1);
mask0=mask;mask0(find(isedge))=0;
holder = sum(logical(mask0(noedge) == 0));
k=1;
if display,
    I{1}=abs(u0).*255./(max(abs(u0(:))));
end
while encore
    un=imfilter(mask0,h,'symmetric');
    f=imfilter(mask0.*u0,h,'symmetric');
    n=find(mask==0&un~=0&~isedge);
    u0(n)=f(n)./un(n);
    mask0(n)=1;
%     mask0=double(un~=0&~isedge);
    encore=double(min(mask0(noedge))==0);
    
    if sum(logical(mask0(noedge) == 0)) == holder,
        encore = 0;
        mask = ~mask0.*2;
    else
        holder = sum(logical(mask0(noedge) == 0));
    end
    k=k+1;
    if display,
        I{k}=abs(u0).*255./(max(abs(u0(:))));
    end
end
if length(find(mask==2))~=0
    u0=cleanOF(u0,(mask~=2));
end
% u0(find(isedge))=NaN;
if display,
    imview(I)
end

% Magic wand

% u0=magicwand(u0,3);
% u0=magicwand(u0,7);
% u0=magicwand(u0,9);
% u0=magicwand(u0,11);
% u0=magicwand(u0,13);
% u0=magicwand(u0,13);

% u0=wapprox(u0);

function u=magicwand(u,w)
u=medfilt2(real(u),w*[1,1])+i*medfilt2(imag(u),w*[1,1]);
u(1:w,1:w)=u(w,w);
u(end-w+(1:w),1:w)=u(end-w+1,w);
u(1:w,end-w+(1:w))=u(w,end-w+1);
u(end-w+(1:w),end-w+(1:w))=u(end-w+1,end-w+1);

function u0=wapprox(u0)
[M,~]=size(u0);
alpha=3;tau=0;type='o';J=floor(log2(M));
[FFTanalysisfilters,FFTsynthesisfilters]=FFTfractsplinefilters_d(M,alpha,tau,type);
w=FFTwaveletanalysis2D(u0,FFTanalysisfilters,J);
w=threshold(w,0.9);
u0=FFTwaveletsynthesis2D(w,FFTsynthesisfilters,J);