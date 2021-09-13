function [m,x] = nsumk(n,k)
% NSUMK Number and listing of non-negative integer n-tuples summing to k
%    M = NSUMK(N,K) where N and K are positive integers returns M=nchoosek(K+N-1,N-1)
%    This is the number of ordered N-tuples of non-negative integers summing to K
%
%    [M,X] = NSUMK(N,K) produces a matrix X with
%    nchoosek(K+N-1,N-1) rows and n columns. Each row comprises
%    non-negative integers summing to k. The ordering of rows follows the
%    same convention as NCHOOSEK, which is undocumented but with some
%    reliability appears to be lexicographic. The reverse of this presumed ordering
%    is a natural way to list coefficients of polynomials in N variables of degree K.
%    As per nchoosek, this syntax is only practical for situations where N is
%    less than about 15.
% 
%  EXAMPLES:   m = nsumk(5,2)   
%              [~,x] = nsumk(5,2) returns a 15 x 5 matrix x in which rows sum to 2

% Peter Cotton (2021). nsumk (https://www.mathworks.com/matlabcentral/fileexchange/28340-nsumk), MATLAB Central File Exchange.

if isscalar(n) && isscalar(k) && nargout<=1
    m = nchoosek(k+n-1,n-1);
elseif isscalar(n) && isscalar(k) && nargout==2
    m = nchoosek(k+n-1,n-1);
    dividers = [zeros(m,1),nchoosek((1:(k+n-1))',n-1),ones(m,1)*(k+n)];
    x = diff(dividers,1,2)-1;
else
    error('nsumk anticipates scalar k and n');
end
