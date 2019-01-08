% Numerical approximation to Poisson's equation over the square [a,b]x[a,b] with
% Dirichlet boundary conditions.  Uses a second order finite differences
% on a cell-centered mesh with m cells in each direction
% Input:
%     ffun : the RHS of poisson equation (i.e. the Laplacian of u).
%     gfun : the boundary function representing the Dirichlet B.C.
%      a,b : the interval defining the square
%        m : the number of cells in each direction of the mesh.
% Ouput:
%        u : the numerical solution of Poisson equation at the mesh points.
%      x,y : the uniform mesh.
%
function [u,x,y] = fd2poissonsp(ffun,gfun,a,b,mx,my)

hx = (b-a)/mx;   % Mesh spacing
hy = (b-a)/my;

[x,y] = meshgrid(a+hx/2:hx:b-hx/2,a+hy/2:hy:b-hy/2); % Uniform mesh

% Compute boundary terms, south, north, east, west
ubs = feval(gfun,x(1,:),y(1,:)-hy/2);  
ubn = feval(gfun,x(my,:),y(my,:)+hy/2); 
ube = feval(gfun,x(:,mx)+hx/2,y(:,mx));   
ubw = feval(gfun,x(:,1)-hx/2,y(:,1)); 

% Evaluate the RHS of Poisson's equation at the interior points.
f = feval(ffun,x,y);

% Adjust f for boundary terms
f(:,1) = f(:,1) - 2*ubw/hx^2;      % West
f(:,mx) = f(:,mx) - 2*ube/hx^2;      % East
f(1,:) = f(1,:) - 2*ubs/hy^2;    % South
f(my,:) = f(my,:) - 2*ubn/hy^2;    % North
f = reshape(f,mx*my,1);

% Create the D2x and D2y matrices

% Sparse matrix version.
D2 = toeplitz(spdiags([-2 1],[0 1],1,mx)); D2(1,1) = -3; D2(mx,mx) = -3;
D2x = 1/hx^2*kron(D2,speye(my));
D2 = toeplitz(spdiags([-2 1],[0 1],1,my)); D2(1,1) = -3; D2(my,my) = -3;
D2y = 1/hy^2*kron(speye(mx),D2);

% Solve the system
u = (D2x + D2y)\f;

% Convert u from a column vector to a matrix to make it easier to work with
% for plotting.
u = reshape(u,my,mx);
 
end



