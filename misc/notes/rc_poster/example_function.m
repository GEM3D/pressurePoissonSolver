u = @(x,y) exp(-40*((x-0.4).^2 + (y-0.4).^2)).*(exp(cos(10*pi*x))-exp(cos(11*pi*y))) + ...
           exp(-10*((x-0.5).^2 + (y-0.55).^2)).*(exp(cos(11*pi*x))-exp(cos(10*pi*y)));
u = chebfun2(u,[0 1 0 1]);
f = laplacian(u);
surf(u), axis equal, view(2)
xlabel('x'), ylabel('y')
set(gca,'YTick',0:0.2:1)
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',14)
colorbar