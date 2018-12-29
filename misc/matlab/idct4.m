function coeffs = idct4(x,dim)
%  IDCT4   Inverse Discrete Cosine Transform Type IV computed using the fast Fourier Transform.
%     X = idct4(x) computes the Inverse Discrete Cosine Transform Type IV (DCT-IV) of the columns of X.
%
%     X = idct4(x,dim) computes the inverse DCT along the dimension specified.
%     if dim = 1 (default) then the inverse DCT is along the columns.
%     if dim = 2 then the inverse DCT is along the rows.
%
%  See also dct4, dct2, idct2, dct, idct, dst, dst2, idst, idst2.

[m,n] = size(x);

if nargin == 1
    dim = 1;
end

% Trivial case (constant):
if ( m <= 1 )
    coeffs = x;
    return
end

% IDCT-IV is a scaled DCT-IV:
if dim == 2
    coeffs = 2/n*dct4(x,dim);
else
    coeffs = 2/m*dct4(x,dim);
end

end