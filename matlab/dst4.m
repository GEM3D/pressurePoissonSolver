function x = dst4(x,dim)
%  DST   Discrete Sine Transform Type IV computed using the fast Fourier Transform.
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

%
% Compute using the DCT-IV.
%
if dim == 1
    x = dct4(x(m:-1:1,:),dim);
    x(2:2:m,:) = -x(2:2:n,:);
elseif dim == 2
    x = dct4(x(:,n:-1:1),dim);
    x(:,2:2:n) = -x(:,2:2:n);
else
    error('dst4:dimUnknown','DST-IV dimension not available, select 1 or 2');
end

end

