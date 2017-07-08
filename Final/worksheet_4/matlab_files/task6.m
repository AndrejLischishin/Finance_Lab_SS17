fid = fopen('../output/task6.txt', 'r');
format_spec = '%f %f %f %f ';
sizeA = [4 Inf];
A = fscanf(fid, format_spec, sizeA);
fclose(fid);
A=A';
n = A(:,1);
x = A(:,2);
y = A(:,3);
z = A(:,4);


plot(log10(n), x,'bo--','Linewidth',2);

hold on

plot(log10(n), y,'r^-','Linewidth',2);

hold on

plot(log10(n), z,':ms','Linewidth',2);





title('task 6 M=64');
ylabel('value');
xlabel('N');
lngd = legend('QMC BB','MC BB','reference value');
set(lngd,'fontsize',4,'Position',[0.6,0.4,0.3,0.1]);
