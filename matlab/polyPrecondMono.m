function y = polyPrecondMono(b,A,coeffs)
% Computes the matrix polynomial product p_n(A)*b where 
% the polynomial coefficients are given in the monomial basis:
%       y = sum_{k=0}^n coeffs(k)*x^k*b 
% The code uses Horner's scheme for evaluating the polynomial times vector
% product.

nc = numel(coeffs);
y = coeffs(nc)*b;
for i = nc-1:-1:1
    y = coeffs(i)*b + A*y;
end

end
