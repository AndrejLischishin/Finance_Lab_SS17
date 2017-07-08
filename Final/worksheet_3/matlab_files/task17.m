fid = fopen('../output/testfunction_task17_error.txt', 'r');
format_spec = '%f %f %f %f %f';
sizeA = [5 Inf];
A = fscanf(fid, format_spec, sizeA);
fclose(fid);
A=A';
n = A(:,1);
x = A(:,2);
y = A(:,3);
z = A(:,4);
u = A(:,5);

loglog(n, x,'bo--','Linewidth',2);

hold on

loglog(n, y,'r^-','Linewidth',2);

hold on

loglog(n, z,':ms','Linewidth',2);

hold on

loglog(n, u,'-g*','Linewidth',2);



title('task 17 errors M=64');
ylabel('error(N)');
xlabel('N');
lngd = legend('QMC RW','QMC BB','MC RW','MC BB');
set(lngd,'fontsize',4,'Position',[0.6,0.8,0.3,0.1]);
