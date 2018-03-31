function y = polyPrecondCheb(b,A,coeffs,intvl)
% Computes the matrix polynomial product p_n(A)*b where 
% the polynomial coefficients, in the Chebyshev basis, are given 
% by coeffs and intvl specifies the interval of the approximation
% If intvl = [-1,1], then what is computed is 
%       y = sum_{k=0}^n coeffs(k)*T_k(A)*b 
% where T_k is the kth degree Chebyshev polynomial. If the interval is more
% general then a simple transformation of A is necessary first.
% The code uses Clenshaw's scheme for evaluating the polynomial time vector
% product.

n = size(A,1);
nc = numel(coeffs);

% The code below assumes the interval is [0 d]. It will not work otherwise
d = intvl(2);
if intvl(1) ~= 0 || d <= 0
    error('Interval must be [0,d], for some d>0')
end

% Clenshaw scheme
bk1 = 0*b; 
bk2 = bk1;
A = 2*((2/d)*A - speye(n));
for k = nc:-1:2
    bk = coeffs(k)*b + A*bk1 - bk2;
    bk2 = bk1; 
    bk1 = bk;
end
y = coeffs(1)*b + 0.5*A*bk1 - bk2;

% Another version of the code that avoids setting intermidiate variables
% every time through the loop, which could be better for parallelization.
% bk1 = 0*b; 
% bk2 = bk1;
% A = 2*((2/d)*A - speye(n));
% for k = nc:-2:3
%     bk2 = coeffs(k)*b + A*bk1 - bk2;
%     bk1 = coeffs(k-1)*b + A*bk2 - bk1;
% end
% if ( mod(nc, 2) )
%     [bk1, bk2] = deal(coeffs(2)*b + A*bk1 - bk2, bk1);
% end
% y = coeffs(1)*b + 0.5*A*bk1 - bk2;

end
