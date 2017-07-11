fid = fopen('../output/task2_error.txt', 'r');
format_spec = '%f %f %f %f %f ';
sizeA = [5 Inf];
A = fscanf(fid, format_spec, sizeA);
fclose(fid);
A=A';
n = A(:,1);
x = A(:,2);
y = A(:,3);
z = A(:,4);
v = A(:,5);


loglog(n, x,'bo--','Linewidth',2);

hold on

loglog(n, y,'--ko','Linewidth',2);

hold on

loglog(n, z,':ms','Linewidth',2);

hold on

loglog(n, v,'r^-','Linewidth',2);


title('M = 64, barrier down-out call, QMC, MC');
ylabel('error[|ref.val. - calc.res.|]');
xlabel('N');
lngd = legend('QMC RW','QMC BB','MC RW','MC BB');
set(lngd,'fontsize',6,'Position',[0.6,0.69,0.3,0.2]);