function x = dst2mm(x,dim)
%  DST   Discrete Sine Transform Type II computed using matrix
%  multiplication.
%     X = dst2(x) computes the Discrete Sine Transform Type II (DST-II) of the columns of X.
%
%     X = dst2(x,dim) computes the DST-II along the dimension specified.
%     if dim = 1 (default) then the DST-II is along the columns.
%     if dim = 2 then the DST is along the rows.
%  
%     See also idst2, dst, idst

if nargin == 1
    dim = 1;
end

[m,n] = size(x);

if dim == 2
    k=1:n;
    A = sin(pi/n*((k-0.5)'*k));
    x = x*A;
elseif dim == 1
    k=1:m;
    A = sin(pi/m*(k'*(k-0.5)));
    x = A*x;
else
    error('dst2mm:dimUnknown','DST-II dimension not available, select 1 or 2');
end

end

