function x = idst2(x,dim)
%  DST   Inverse Discrete Sine Transform Type II computed using the fast Fourier Transform.
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

%
% Compute using the inverse DCT-II.
%

if dim == 1
    x(1:2:m,:) = -x(1:2:m,:);
    x = idct2(x(m:-1:1,:),1);
    x(1:2:m,:) = -x(1:2:m,:);
    x = x(m:-1:1,:);
elseif dim == 2
    x(:,1:2:n) = -x(:,1:2:n);
    x = idct2(x(:,n:-1:1),2);
    x(:,1:2:n) = -x(:,1:2:n);
    x = x(:,n:-1:1);
else
    error('idst2:dimUnknown','Inverse DST-II dimension not available, select 1 or 2');
end

end

