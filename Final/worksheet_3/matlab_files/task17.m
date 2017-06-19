fid = fopen('../output/testfunction_task17_error.txt', 'r');
format_spec = '%f %f %f %f %f';
sizeA = [5 Inf];
A = fscanf(fid, format_spec, sizeA);
fclose(fid);
A=A'
n = A(:,1);
x = A(:,2);
y = A(:,3);
z = A(:,4);
u = A(:,5);

loglog(n, x);

hold on

loglog(n, y);

hold on

loglog(n, z);

hold on

loglog(n, u);


title('task 17 errors M=64');
ylabel('error(N)');
xlabel('N');
