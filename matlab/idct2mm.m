function coeffs = idct2mm(x,dim)
%  IDCT2   Inverse Discrete Cosine Transform Type II computed using
%  matrix multiplication.
%     X = idct2mm(x) computes the Inverse Discrete Cosine Transform Type II (DCTII) of the columns of X.
%
%     X = idct2mm(x,dim) computes the inverse DCT along the dimension specified.
%     if dim = 1 (default) then the inverse DCT is along the columns.
%     if dim = 2 then the inverse DCT is along the rows.
%
%  See also dct4, idct2, dct, idct, dst, dst2, idst, idst2.

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
    A = cos(pi/n*(k'*(k+0.5)));
    A(1,:) = 0.5;
    coeffs = 2/n*(x*A);
else
    k=0:m-1;
    A = cos(pi/m*((k+0.5)'*k));
    A(:,1) = 0.5;
    coeffs = 2/m*(A*x);
end

end
