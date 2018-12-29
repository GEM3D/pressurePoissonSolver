function coeffs = dct2mm(x,dim)
%  DCT2   Discrete Cosine Transform Type II computed using matrix
%  multiplication
%     X = dct2(x) computes the Discrete Cosine Transform Type II (DCT-II) of the columns of X.
%
%     X = dct2(x,dim) computes the DCT along the dimension specified.
%     if dim = 1 (default) then the DCT is along the columns.
%     if dim = 2 then the DCT is along the rows.
%
%  See also idct2, dct, idct, dst, dst2, idst, idst2.

[m,n] = size(x);

if nargin == 1
    dim = 1;
end

% Trivial case (constant):
if ( m <= 1 )
    coeffs = x;
    return
end

if dim == 2
    k=0:n-1;
    A = cos(pi/n*((k+0.5)'*k));
    coeffs = x*A;
else
    k=0:m-1;
    A = cos(pi/m*(k'*(k+0.5)));
    coeffs = A*x;
end

end