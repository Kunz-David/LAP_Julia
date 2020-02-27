function J=ShiftedLinear_Interp_2D(ux,uy,I)
% USAGE    : J=ShiftedLinear_Interp_2D(ux,uy,I)
% FUNCTION : Interpolates the 2D image I at the (in general, noninteger) 
% positions x - ux and y - uy. ux and uy are 2D matrices with 
% same size and represent displacements in x and y, respectively. The 
% function combines shifted-linear interpolation with matlab's inbuilt 
% interp2 (with Bilinear input). The result is a fast interpolation that 
% is more accurate than standard Trilinear interpolation (approximately 
% equivalent to Key's cubic interpolation). For more details see:
%
% T. Blu, P. Thevenaz, and M. Unser, "Linear interpolation revitalized,"
%   IEEE Trans. Image Processing, vol. 13, no. 5, pp. 710–719, 2004.
%
% DATE     : 23 December 2015
% AUTHOR   : Christopher Gilliam, email: dr.christopher.gilliam@ieee.org
%
% Note: Based on the more general imshift/interp function written by 
% Thierry Blu, mailto:thierry.blu@m4x.org
%

[M,N]=size(I);

% Generate x, y and z:
x=(1:M)'*ones(1,N);
y=ones(M,1)*(1:N);

% Parameters for shifted-linear interpolation:
tau=1/2*(1-sqrt(3)/3);
L1=floor(-1+tau);      
L2=ceil(1+tau);

x1 = x - ux;
y1 = y - uy;

% Minimum and maximum row index needed in the interpolation formula
k0=floor(min(x1(:))-L2+1);
k1=floor(max(x1(:))-L1);
l0=floor(min(y1(:))-L2+1);
l1=floor(max(y1(:))-L1);

% Smallest box enclosing the image and the (x,y) positions
kk0=min(k0,1);
kk1=max(k1,M);
ll0=min(l0,1);
ll1=max(l1,N);

% For 3D images:
% Image extension to fill the unknown pixels
exttype='symmetric';                            % options are: 
                                                % 'circular' (Pad with circular repetition of elements within the dimension), 
                                                % 'replicate' (Pad by repeating border elements of array), 
                                                % 'symmetric' (Pad array with mirror reflections of itself), 

I0=ext(I,exttype,[1-kk0 kk1-M 1-ll0 ll1-N]);
I0=I0(1-kk0+(k0:k1),1-ll0+(l0:l1));
[a0,b0]=size(I0);

x = repmat((k0:k1)',[1,b0]);
y = repmat((l0:l1), [a0,1]);

% Prefiltering image I0 so as to have shifted-linear interpolation:
z0=tau/(1-tau);
% along columns first
I0=1/(1-tau)*filtering(1,[1 z0],I0,'causal');
% then lines
I0=I0.';
I0=1/(1-tau)*filtering(1,[1 z0],I0,'causal');
I0=I0.';

% Matlab's inbuilt interpolation command:
J = interp2(y,x,I0,y1-tau,x1-tau,'linear');

%*************************************************
%*************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Embedded functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J=ext(I,exttype,extsize)
[a,b]=size(I);
newa=a+extsize(1)+extsize(2);
newb=b+extsize(3)+extsize(4);

% Extend in x:
if extsize(1)>extsize(2)
%     J=wextend(2,exttype,I,extsize(1),'bn');
    J=padarray(I,[extsize(1),0,0],exttype);
    J=J(1:newa,:,:);
else
%     J=wextend(2,exttype,I,extsize(2),'bn');
    J=padarray(I,[extsize(2),0,0],exttype);
    J=J(end+(1:newa)-newa,:,:);
end

% Extend in y:
if extsize(3)>extsize(4)
%     J=wextend(2,exttype,J,extsize(3),'nb');
    J=padarray(J,[0,extsize(3),0],exttype);
    J=J(:,1:newb,:);
else
%     J=wextend(2,exttype,J,extsize(4),'nb');
    J=padarray(J,[0,extsize(4),0],exttype);
    J=J(:,end+(1:newb)-newb,:);
end


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J=filtering(b,a,I,type)
switch type
    case 'causal'
        J=bsxfun(@plus,filter(b,a,bsxfun(@minus,I,I(1,:,:))),I(1,:,:)*sum(b)/sum(a));
    case 'anticausal'
        I=I(end:-1:1,:,:);
        J=bsxfun(@plus,filter(b,a,bsxfun(@minus,I,I(1,:,:))),I(1,:,:)*sum(b)/sum(a));
        J=J(end:-1:1,:,:);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%