function est_noise_variance = estimation_noise_variance (y)
% USAGE    : est_noise_variance = estimation_noise_variance (y)
% FUNCTION : Estimate the noise variance of an image y, by evaluating 
%            the Median of its Absolute Difference
%
% DATE     : 28 July 2015
% AUTHOR   : Thierry Blu, the Chinese University of Hong kong, Shatin, Hong Kong
%            mailto:thierry.blu@m4x.org

[M,N,P] = size(y);

if P == 1,
    dy=abs(imfilter(y,[0 1 0;1 -4 1;0 1 0]/sqrt(20)));
    est_noise_variance=(median(dy(:))/(erfcinv(0.5)*sqrt(2)))^2;
%     dy=diff(diff(y,2,1),2,2);
%     est_noise_variance=(median(abs(dy(:)))/(6*sqrt(2)*erfcinv(0.5)))^2;
else
    dy=diff(diff(diff(y,2,1),2,2),2,3);
    est_noise_variance=(median(abs(dy(:)))/(15*sqrt(2)*erfcinv(0.5)))^2;
end
