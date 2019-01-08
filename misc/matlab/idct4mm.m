function coeffs = idct4mm(x,dim)
%  IDCT4   Inverse Discrete Cosine Transform Type IV computed using matrix
%  multiplication.
%     X = idct4mm(x) computes the Inverse Discrete Cosine Transform Type IV (DCT-IV) of the columns of X.
%
%     X = idct4mm(x,dim) computes the inverse DCT along the dimension specified.
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
    coeffs = 2/n*dct4mm(x,dim);
else
    coeffs = 2/m*dct4mm(x,dim);
end

end