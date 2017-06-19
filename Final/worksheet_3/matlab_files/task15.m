fid = fopen('../output/testfunction_task15_error.txt', 'r');
format_spec = '%f %f %f';
sizeA = [3 Inf];
A = fscanf(fid, format_spec, sizeA);
fclose(fid);
A=A';
n = A(:,1);
x = A(:,2);
y = A(:,3);

loglog(n, x, 'bo--','Linewidth',2);

hold on

loglog(n, y, 'r^-','Linewidth',2);



title('task 15 errors 16d');
ylabel('error(N)');
xlabel('N');


lngd = legend('CC RW','CC BB');
