fid = fopen('../output/testfunction_task13_error_d2.txt', 'r');
format_spec = '%f %f %f %f %f %f %f';
sizeA = [7 Inf];
A = fscanf(fid, format_spec, sizeA);
fclose(fid);
A=A'
n = A(:,1);
x = A(:,2);
y = A(:,3);
z = A(:,4);
u = A(:,5);
v = A(:,6);
w = A(:,7);

loglog(n, x, "0;QMC;");

hold on

loglog(n, y, "5;MC;");

hold on

loglog(n, z, "4;FG trap rule;");

hold on

loglog(n, u, "3;FG CC;");

hold on

loglog(n, v, "2;SG trap rule;");

hold on

loglog(n, w, "1;SG CC;");

title('task 13 errors 2d');
ylabel('error(N)');
xlabel('N');