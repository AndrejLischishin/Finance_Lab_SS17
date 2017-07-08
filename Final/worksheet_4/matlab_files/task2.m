fid = fopen('../output/task2.txt', 'r');
format_spec = '%f %f %f %f %f %f ';
sizeA = [6 Inf];
A = fscanf(fid, format_spec, sizeA);
fclose(fid);
A=A';
n = A(:,1);
x = A(:,2);
y = A(:,3);
z = A(:,4);
v = A(:,5);
w = A(:,6);


plot(log10(n), x,'bo--','Linewidth',2);

hold on

plot(log10(n), y,'--ko','Linewidth',2);

hold on

plot(log10(n), z,':ms','Linewidth',2);

hold on

plot(log10(n), v,'r^-','Linewidth',2);

hold on

plot(log10(n), w,'--mo','Linewidth',2);


title('task 2 M=64');
ylabel('value');
xlabel('N');
lngd = legend('QMC RW','QMC BB','MC RW','MC BB','reference value');
set(lngd,'fontsize',4,'Position',[0.6,0.4,0.3,0.1]);