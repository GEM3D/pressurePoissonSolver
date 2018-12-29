function x = dst2(x,dim)
%  DST   Discrete Sine Transform Type II computed using the fast Fourier Transform.
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

%
% Compute using the DCT-II.
%

if dim == 1
    x(1:2:m,:) = -x(1:2:m,:);
    x = dct2(x(m:-1:1,:),1);
    x(1:2:m,:) = -x(1:2:m,:);
    x = x(m:-1:1,:);
elseif dim == 2
    x(:,1:2:n) = -x(:,1:2:n);
    x = dct2(x(:,n:-1:1),2);
    x(:,1:2:n) = -x(:,1:2:n);
    x = x(:,n:-1:1);
else
    error('dst2:dimUnknown','DST-II dimension not available, select 1 or 2');
end

end

