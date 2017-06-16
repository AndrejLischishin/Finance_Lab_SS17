fid = fopen('../output/number_of_points_SG_FG.txt', 'r');
format_spec = '%f %f %f';
sizeA = [3 Inf];
A = fscanf(fid, format_spec, sizeA);
fclose(fid);
A=A'
d = A(:,1);
x = A(:,2);
y = A(:,3);

semilogy(d, x, "2;Sparse Grid;");

hold on

semilogy(d, y, "3;Full Grid;");

title('Number of points in Full and Sparse Grid of Level l=4');
ylabel('Number of points');
xlabel('d');