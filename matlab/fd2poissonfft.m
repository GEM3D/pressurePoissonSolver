% Numerical approximation to Poisson's equation over the square [a,b]x[a,b] with
% Dirichlet boundary conditions.  Uses a second order finite differences
% on a cell-centered mesh with mx-by-my cells in the directions.  Solves
% the equations using the FFT (specifically the fast sine transform II -
% see http://en.wikipedia.org/wiki/Discrete_sine_transform).
% Input:
%     ffun : the RHS of poisson equation (i.e. the Laplacian of u).
%     gfun : the boundary function representing the Dirichlet B.C.
%      a,b : the interval defining the square
%       mx : the number of cells in x direction of the mesh.
%       my : the number of cells in y direction of the mesh.
% Ouput:
%        u : the numerical solution of Poisson equation at the mesh points.
%      x,y : the uniform mesh.
%
function [u,x,y] = fd2poissonfft(ffun,gfun,a,b,mx,my)

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

jx = (1:mx);
jy = (1:my)';
fhat = dst2(dst2(f,1),2);
den = bsxfun(@plus, -4/hx^2*sin(jx*pi/(2*mx)).^2,-4/hy^2*sin(jy*pi/(2*my)).^2);

uhat = fhat./den;
u = idst2(idst2(uhat,1),2);
 
end



