function basis=loadbasis_asymm(K)
% function generates the elements for the filter basis. Case 2 is the
% Gaussian filter basis defined in the ICASSP2015 and ICIP2015 papers. The
% parameter K defined the half-size of the filtes in x and y (i.e. the
% filters are 2*K(1) +1 by 2*K(2) +1 in size)
%

% 
btype=2;

switch btype
    case 1
        % old filter definitions (do not use)
        basis=zeros(5,5,3);
        basis(:,:,1)=[0 0 1 0 0;0 15 30 15 0;1 30 0 30 1;0 15 30 15 0;0 0 1 0 0];
        basis(:,:,2)=[0 0 2 0 0;0 15 30 11 0;0 4 0 -4 0;0 -11 -30 -15 0;0 0 -2 0 0];
        basis(:,:,3)=basis(:,:,2)';
        basis=reshape(basis,25,3);
    case 2
        if length(K) == 1,
            K(1) = K;
            K(2) = K;
        end
        s1=(K(1)+2)/4;
        s2=((K(2)+2)/4);
        g=@(x,s)exp(-x.*x/2/s^2);
        
        if length(K) > 1,
            K1=ceil(K(1));
            K2=ceil(K(2));
        else
            K1=ceil(K);
            K2 = K1;
        end
        k=-K1:K1;k=k'*ones(1,2*K2+1);
        l=ones(2*K1+1,1)*(-K2:K2);
        
        basis=zeros(2*K1+1,2*K2+1,6);
        basis(:,:,1)=g(k,s1).*g(l,s2);basis(:,:,1)=basis(:,:,1)/sum(sum(basis(:,:,1)));
        basis(:,:,2)=g(k,s1).*l.*g(l,s2);basis(:,:,2)=basis(:,:,2)/sum(sum(l.*basis(:,:,2)));
        basis(:,:,3)=g(k,s1).*k.*g(l,s2);basis(:,:,3)=basis(:,:,3)/sum(sum(k.*basis(:,:,3)));  
        basis(:,:,4)=g(k,s1).*(k.*k+l.*l).*g(l,s2)/(s1*s2);
        basis(:,:,4)=basis(:,:,4)-basis(:,:,1)*sum(sum(basis(:,:,4)))/sum(sum(basis(:,:,1)));
        basis(:,:,5)=g(k,s1).*(k.*k-l.*l).*g(l,s2)/(s1*s2);
        basis(:,:,6)=g(k,s1).*(k.*l).*g(l,s2)/(s1*s2);
        basis=reshape(basis,(2*K1+1)*(2*K2+1),6);
        
    case 3
        % generate a single Gaussian (do not use)
        s1=1+(K(1)-2)/4;
        s2=1+(K(2)-2)/4;
        g=@(x,s)exp(-x.*x/2/s^2);
        
        if length(K) > 1,
            K1=ceil(K(1));
            K2=ceil(K(2));
        else
            K1=ceil(K);
            K2 = K1;
        end
        k=-K1:K1;k=k'*ones(1,2*K2+1);
        l=ones(2*K1+1,1)*(-K2:K2);
        basis(:,:,1)=g(k,s1).*g(l,s2);basis(:,:,1)=basis(:,:,1)/sum(sum(basis(:,:,1)));
        basis=reshape(basis,(2*K1+1)*(2*K2+1),1);
        
end
end