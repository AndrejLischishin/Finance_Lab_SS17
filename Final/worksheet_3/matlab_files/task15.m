fid = fopen('../output/testfunction_task15_error.txt', 'r');
format_spec = '%f %f %f';
sizeA = [3 Inf];
A = fscanf(fid, format_spec, sizeA);
fclose(fid);
A=A'
n = A(:,1);
x = A(:,2);
y = A(:,3);

loglog(n, x, "0;CC Random Walk;");

hold on

loglog(n, y, "5;CC Brownian Bridge;");



title('task 15 errors 16d');
ylabel('error(N)');
xlabel('N');
