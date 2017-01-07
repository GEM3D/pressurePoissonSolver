function coeffs = idct2(x,dim)
%  IDCT2   Inverse Discrete Cosine Transform Type II computed using the fast Fourier Transform.
%     X = idct2(x) computes the Inverse Discrete Cosine Transform Type II (DCTII) of the columns of X.
%
%     X = idct2(x,dim) computes the inverse DCT along the dimension specified.
%     if dim = 1 (default) then the inverse DCT is along the columns.
%     if dim = 2 then the inverse DCT is along the rows.
%
%  See also dct2.

% % IDCT-II is a scaled DCT-III.
% y = ( 2 / n ) * chebfun.dct(u, 3); 
% 
% y = chebtech1.coeffs2vals( u );    
% y = y(end:-1:1,:); 

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
    % Account for constant term
    x(:,1) = .5*x(:,1); 

    j = 1:m;
    k = n:-1:1;
    % Mirror the values for FFT:
    tmp = [x ones(m,1) x(:,end:-1:2)];
    scale = n;

    % Pre-compute the weight vector:
    w = (exp(-1i*(0:2*n-1)*pi/(2*n))/2);
    w([1, n+1]) = [2*w(1) 0];
    w(n+2:end) = -w(n+2:end);

else
    j = m:-1:1;
    k = 1:n;

    % Account for constant term
    x(1,:) = .5*x(1,:); 

    % Mirror the values for FFT:
    tmp = [x ; ones(1, n) ; x(end:-1:2,:)];
    scale = m;
    
    % Pre-compute the weight vector:
    w = (exp(-1i*(0:2*m-1)*pi/(2*m))/2).';
    w([1, m+1]) = [2*w(1); 0];
    w(m+2:end) = -w(m+2:end);
end

% Apply the weight vector:
tmp = bsxfun(@times, tmp, w);

coeffs = fft(tmp,[],dim);

% Truncate, flip the order, and scale:
% coeffs = (2/scale)*coeffs(j,k);
coeffs = (2/scale)*coeffs(1:m,1:n);

% Post-process:
if ( isreal(x) )  
    % Real-valued case:
    coeffs = real(coeffs);
elseif ( isreal(1i*x) )  
    % Imaginary-valued case:
    coeffs = 1i*imag(coeffs);
end

end

% 
% % Scale the coefficient for the constant term:
% coeffs = (scale)/2*coeffs;
% 
% % Post-process:
% if ( isreal(x) )  
%     % Real-valued case:
%     coeffs = real(coeffs);
% elseif ( isreal(1i*x) )  
%     % Imaginary-valued case:
%     coeffs = 1i*imag(coeffs);
% end
