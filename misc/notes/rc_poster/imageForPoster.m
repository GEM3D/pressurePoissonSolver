%% Plot 3-by-3 grid
[xx,yy] = meshgrid(linspace(0,1,4));
plot(xx,yy,'k-','LineWidth',2); 
hold on;
plot(xx.',yy.','k-','LineWidth',2);
hold off
axis equal
xlim([-0.05 1.05]);
ylim([-0.05 1.05]);
axis off
set(gcf,'Color',[1 1 1])

%% Plot 4-by-4 grid
h = 1/(4*3);
[xxg,yyg] = meshgrid(h/2:h:1-h/2);
plot(xx,yy,'k-','LineWidth',2); 
hold on;
plot(xx.',yy.','k-','LineWidth',2);
plot(xxg,yyg,'b.','MarkerSize',20); 
plot(xxg(1,:),1/3 + 0*xxg(1,:),'ro','MarkerSize',14,'LineWidth',1)
plot(xxg(1,:),2/3 + 0*xxg(1,:),'ro','MarkerSize',14,'LineWidth',1)
hold off
axis equal
xlim([-0.05 1.05]);
ylim([-0.05 1.05]);
axis off
set(gcf,'Color',[1 1 1])

%% Plot 6-by-6 grid
[xx,yy] = meshgrid(linspace(1,2,7),linspace(0,1,7));
h = 1/(4*6);
[xxg,yyg] = meshgrid(1+h/2:h:2-h/2,h/2:h:1-h/2);
hold on;
plot(xx,yy,'k-','LineWidth',2); 
plot(xx.',yy.','k-','LineWidth',2);
plot(xxg,yyg,'b.','MarkerSize',20); 
plot(1+0*yyg(:,1),yyg(:,1),'ro','MarkerSize',14,'LineWidth',1)
plot(1+1/6+0*yyg(:,1),yyg(:,1),'ro','MarkerSize',14,'LineWidth',1)
plot(xxg(1,:),0.5 + 0*xxg(1,:),'ro','MarkerSize',14,'LineWidth',1)
plot(xxg(1,:),0.5+1/6 + 0*xxg(1,:),'ro','MarkerSize',14,'LineWidth',1)
plot(xxg(1,:),0.5-1/6 + 0*xxg(1,:),'ro','MarkerSize',14,'LineWidth',1)
hold off
axis equal
xlim([0.75 1.25]);
ylim([0.3 0.7]);
axis off
set(gcf,'Color',[1 1 1])

