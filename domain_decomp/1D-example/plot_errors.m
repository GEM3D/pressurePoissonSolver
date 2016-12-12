function plot_errors(N,e,m)
% PLOT_ERRORS plots errors
%
% PLOT_ERRORS(N,E,M) plots errors E on N grids cells, over M domains.
% where N is the total number of grid cells.
%
% 

clf;

N = N(:);
e = e(:);

loglog(N,e,'r.-','markersize',30,'linewidth',2);

mb = polyfit(log(N),log(e),1);
hold on;
loglog(N,exp(polyval(mb,log(N))),'k--');


str = sprintf('Slope = %6.1f',mb(1));
h = legend({'Errors',str},'location','southwest');
set(h,'fontsize',16);
xlabel('N','fontsize',16);
ylabel('Error','fontsize',16);
set(gca,'fontsize',16);
tstr = sprintf('Errors over %d domain(s)',m);
title(tstr,'fontsize',18);

convrate = log(e(1:end-1)./e(2:end))/log(2);

fprintf('\n');
fprintf('%4s %12s %12s\n','N','Error','Conv. Rate');
write_line(30);
for i = 1:length(e),
    if (i == 1)        
        fprintf('%4d %12.2e %12s\n',N(i),e(i),'--  ');
    else
        fprintf('%4d %12.2e %12.4f\n',N(i),e(i),convrate(i-1));
    end
end
write_line(30);
        

hold on;

end

function write_line(n)

for i = 1:n,
    fprintf('-');
end
fprintf('\n');

end