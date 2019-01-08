function x = idst2mm(x,dim)
%  DST   Inverse Discrete Sine Transform Type II computed using matrix
%  multiplication.
%     X = idst2(x) computes the inverse Discrete Sine Transform Type II (DST-II) of the columns of X.
%
%     X = idst2(x,dim) computes the inverse DST-II along the dimension specified.
%     if dim = 1 (default) then the inverse DST-II is along the columns.
%     if dim = 2 then the inverse DST is along the rows.
%  
%     See also dst2, dst, idst

if nargin == 1
    dim = 1;
end

[m,n] = size(x);

if dim == 2
    k=0:n-1;
    A = sin(pi/n*((k+1)'*(k+0.5)));
    A(n,:) = 0.5*(-1).^k;
    x = 2/n*(x*A);
elseif dim == 1
    k=0:m-1;
    A = sin(pi/m*((k+0.5)'*(k+1)));
    A(:,m) = 0.5*(-1).^k;
    x = 2/m*(A*x);
else
    error('idst2mm:dimUnknown','Inverse DST-II dimension not available, select 1 or 2');
end

end

