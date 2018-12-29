function coeffs = dct2(x,dim)
%  DCT2   Discrete Cosine Transform Type II computed using the fast Fourier Transform.
%     X = dct2(x) computes the Discrete Cosine Transform Type II (DCT-II) of the columns of X.
%
%     X = dct2(x,dim) computes the DCT along the dimension specified.
%     if dim = 1 (default) then the DCT is along the columns.
%     if dim = 2 then the DCT is along the rows.
%
%  See also idct2, dct, idct, dst, dst2, idst, idst2.

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
    j = (0:n-1);
    % Mirror the values for FFT:
    tmp = [x x(:,n:-1:1)];
    scale = n;
else
    j = (0:m-1)';
    % Mirror the values for FFT:
    tmp = [x; x(m:-1:1,:)];
    scale = m;
end

% Pre-compute the weight vector:
w = 2*exp(1i*j*pi/(2*scale));

% tmp = [x; x(m:-1:1,:)];
coeffs = ifft(tmp,[],dim);

% Truncate, flip the order, and multiply the weight vector:
coeffs = bsxfun(@times, w, coeffs(1:m,1:n));

% Scale the coefficient for the constant term:
coeffs = (scale)/2*coeffs;

% Post-process:
if ( isreal(x) )  
    % Real-valued case:
    coeffs = real(coeffs);
elseif ( isreal(1i*x) )  
    % Imaginary-valued case:
    coeffs = 1i*imag(coeffs);
end

end