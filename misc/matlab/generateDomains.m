%%
n = 4;
m = 4;
hg = 1/(m-1);
[xd,yd] = meshgrid(linspace(0,1,m));
h = hg/n;
t = linspace(0,hg,n+1); t = t(1:end-1) + 0.5*h;
[xg,yg] = meshgrid(t);

%%
plotDomains = 1;
plotGrid = 1;
plotGamma = 1;
labelDomains = 0;


%%
if plotDomains == 1
    plot(xd,yd,'k-',xd',yd','k-','LineWidth',2)
    hold on
end

if labelDomains == 1
    count = 1;
    for j=1:m-1
        for k=1:m-1
            text(0.5*(xd(j,k)+xd(j+1,k+1)),0.5*(yd(j+1,k+1)+yd(j,k)),sprintf('$\\Omega_{%d}$',count),'Interpreter','Latex','FontSize',16)
            count = count + 1;
        end
    end
end

if plotGrid == 1
    for i = 1:m-1
        for j = 1:m-1
            plot(xg+xd(i,j),yg+yd(i,j),'b.','MarkerSize',10)
            hold on;
        end
    end
end

if plotGamma == 1
    for i = 1:m-1
        for j = 1:m-2
            plot(xg(:,end)+0.5*h+xd(i,j),yg(:,end)+yd(i,j),'ro','MarkerSize',6)
            hold on;
        end
    end
    for i = 1:m-2
        for j = 1:m-1
            plot(xg(end,:)+xd(i,j),yg(end,:)+0.5*h+yd(i,j),'ro','MarkerSize',6)
            hold on;
        end
    end
end

axis off
axis equal
hold off
set(gcf,'Position',[733   765   507   333]);
shg

