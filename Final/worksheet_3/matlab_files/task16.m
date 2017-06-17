fid = fopen('../output/testfunction_task16_error.txt', 'r');
format_spec = '%f %f %f %f %f %f %f %f %f %f %f %f %f';
sizeA = [13 Inf];
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
m = A(:,8);
l = A(:,9);
k = A(:,10);
i = A(:,11);
j = A(:,12);
h = A(:,13);

loglog(n, x, "QMC RW;");

hold on

loglog(n, y, "QMC BB;");

hold on

loglog(n, z, "MC RW;");

hold on

loglog(n, u, "MC BB;");

hold on

loglog(n, v, "FG TRAPEZ. RW;");

hold on

loglog(n, w, "FG TRAPEZ. BB;");

hold on

loglog(n, m, "FG CC RW;");

hold on

loglog(n, l, "FG CC BB");

hold on

loglog(n, k, "SG TRAPEZ RW;");

hold on

loglog(n, i, "SG TRAPEZ BB;");

hold on

loglog(n, j, "SG CC RW;");

hold on

loglog(n, h, "SG CC BB;");

title('task 16 errors');
ylabel('error(N)');
xlabel('N');
