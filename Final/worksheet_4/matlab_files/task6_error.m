fid = fopen('../output/task6_error.txt', 'r');
format_spec = '%f %f %f ';
sizeA = [3 Inf];
A = fscanf(fid, format_spec, sizeA);
fclose(fid);
A=A';
n = A(:,1);
x = A(:,2);
y = A(:,3);


loglog(n, x,'bo--','Linewidth',2);

hold on

loglog(n, y,'r^-','Linewidth',2);


title('M = 64, lookback, QMC, MC');
ylabel('error[|ref.val. - calc.val.|]');
xlabel('N');
lngd = legend('QMC BB','MC BB');
set(lngd,'fontsize',6,'Position',[0.6,0.69,0.3,0.2]);
