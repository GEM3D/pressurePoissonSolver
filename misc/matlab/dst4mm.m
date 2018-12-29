function x = dst4mm(x,dim)
%  DST   Discrete Sine Transform Type IV computed using matrix
%  multiplication.
%     X = dst4(x) computes the Discrete Sine Transform Type IV (DST-IV) of the columns of X.
%
%     X = dst4(x,dim) computes the DST-IV along the dimension specified.
%     if dim = 1 (default) then the DST-IV is along the columns.
%     if dim = 2 then the DST-IV is along the rows.
%  
%  See also dct4, idct4, dct2, idct2, dct, idct, dst, dst2, idst, idst2,
%  idst4.

if nargin == 1
    dim = 1;
end

[m,n] = size(x);

if dim == 2
    k=0:n-1;
    A = sin(pi/n*((k+0.5)'*(k+0.5)));
    x = x*A;
elseif dim == 1
    k=0:m-1;
    A = sin(pi/m*((k+0.5)'*(k+0.5)));
    x = A*x;
else
    error('dst4:dimUnknown','DST-IV dimension not available, select 1 or 2');
end

end

