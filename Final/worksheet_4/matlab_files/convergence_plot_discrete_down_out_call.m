fid = fopen('../output/convergence_plot_discrete_down_out_call.txt', 'r');
format_spec = '%f %f %f %f %f';
sizeA = [5 Inf];
A = fscanf(fid, format_spec, sizeA);
fclose(fid);
d = A(1,:);
x = A(2,:);
y = A(3,:);
z = A(4,:);
w = A(5,:);

loglog(d, x, "2;M=4;");

hold on

loglog(d, y, "3;M=64;");

hold on

loglog(d, z, "4;M=256;");

hold on

loglog(d, w, "5;M=1024;");

title('Absolute error for discrete Down-out call');
ylabel('Error');
xlabel('n');