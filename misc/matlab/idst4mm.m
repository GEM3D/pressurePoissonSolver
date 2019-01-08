function coeffs = idst4mm(x,dim)
%  IDST4   Inverse Discrete Sine Transform Type IV computed using matrix
%  multiplication.
%     X = idst4mm(x) computes the Inverse Discrete Sine Transform Type IV (DST-IV) of the columns of X.
%
%     X = idst4mm(x,dim) computes the inverse DST along the dimension specified.
%     if dim = 1 (default) then the inverse DST is along the columns.
%     if dim = 2 then the inverse DST is along the rows.
%
%  See also dct4, idct4, dct2, idct2, dct, idct, dst, dst2, idst, idst2,
%  dst4.

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
if dim == 1
    coeffs = 2/m*dst4mm(x,dim);
elseif dim == 2
    coeffs = 2/n*dst4mm(x,dim);
else
    error('idst4mm:dimUnknown','IDST-IV dimension not available, select 1 or 2');
end

end