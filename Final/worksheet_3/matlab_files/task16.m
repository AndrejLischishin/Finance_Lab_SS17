fid = fopen('../output/testfunction_task16_error.txt', 'r');
format_spec = '%f %f %f %f %f %f %f %f %f %f %f %f %f';
sizeA = [13 Inf];
A = fscanf(fid, format_spec, sizeA);
fclose(fid);
A=A';
n = A(:,1);
x = A(:,2);
y = A(:,3);
z = A(:,4);
u = A(:,5);
v = A(:,6);
w = A(:,7);
m = A(:,8);
l = A(:,9);
k = A(:,10);
i = A(:,11);
j = A(:,12);
h = A(:,13);

subplot(2,1,1);
loglog(n, x,'bo--','Linewidth',2);

hold on

loglog(n, z,'r^-','Linewidth',2);

hold on

loglog(n, v,':ms','Linewidth',2);

hold on

loglog(n, m,'-g*','Linewidth',2);

hold on

loglog(n, k,'--mo','Linewidth',2);

hold on

loglog(n, j,'--ko','Linewidth',2);
title('All methods using RW')
ylabel('error(N)');
xlabel('N');

lngd = legend('QMC','MC','FG TRAP','FG CC','SG TRAP','SG CC');
set(lngd,'fontsize',7,'Position',[0.6,0.73,0.3,0.1]);


%%%%%%%%%%%%%%%%%%%%%%
%%%%%%next plot%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%

subplot(2,1,2);
loglog(n, y,'bo--','Linewidth',2);

hold on

loglog(n, u,'r^-','Linewidth',2);

hold on

loglog(n, w,':ms','Linewidth',2);

hold on

loglog(n, l,'-g*','Linewidth',2);

hold on

loglog(n, i,'--mo','Linewidth',2);

hold on

loglog(n, h,'--ko','Linewidth',2);
title('All methods using BB')
ylabel('error(N)');
xlabel('N');
lngd = legend('QMC','MC','FG TRAP','FG CC','SG TRAP','SG CC');
set(lngd,'fontsize',4,'Position',[0.6,0.27,0.3,0.1]);
