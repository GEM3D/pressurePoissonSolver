function coeffs = dct4(x,dim)
%  DCT4   Discrete Cosine Transform Type IV computed using the fast Fourier Transform.
%     X = dct4(x) computes the Discrete Cosine Transform Type IV (DCT-IV) of the columns of X.
%
%     X = dct4(x,dim) computes the DCT along the dimension specified.
%     if dim = 1 (default) then the DCT is along the columns.
%     if dim = 2 then the DCT is along the rows.
%
%  See also idct4, dct2, idct2, dct, idct, dst, dst2, idst, idst2.

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
    v = zeros(m,2*n);
    v(:,2:2:2*n) = x;
    n = 2*n;
    % Mirror the values for FFT:
    tmp = [v, ones(m,1), v(:,n:-1:2)];
    % Pre-compute the weight vector:
    j = (0:2*n-1);
    w = exp(-1i*j*pi/(2*n))/2;
    w([1, n+1]) = [2*w(1) 0];
    w(1,n+2:end) = -w(1,n+2:end);
    % Apply the weight vector:
    coeffs = fft(bsxfun(@times, tmp, w),[],dim);
    % Extract the coefficients
    coeffs = coeffs(:,1:n/2);
else
    v = zeros(2*m,n);
    v(2:2:2*m,:) = x;
    m = 2*m;
    % Mirror the values for FFT:
    tmp = [v; ones(1,n); v(m:-1:2,:)];
    % Pre-compute the weight vector:
    j = (0:2*m-1)';
    w = exp(-1i*j*pi/(2*m))/2;
    w([1, m+1]) = [2*w(1); 0];
    w(m+2:end) = -w(m+2:end);
    % Apply the weight vector:
    coeffs = fft(bsxfun(@times, tmp, w),[],dim);
    % Truncate and flip the order:
    coeffs = coeffs(1:m/2,:);
end

% Post-process:
if ( isreal(x) )  
    % Real-valued case:
    coeffs = real(coeffs);
elseif ( isreal(1i*x) )  
    % Imaginary-valued case:
    coeffs = 1i*imag(coeffs);
end

end