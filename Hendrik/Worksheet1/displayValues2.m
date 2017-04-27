fid = fopen('values2.txt', 'r');
data = fscanf(fid, '%f');
fclose(fid);
[vals, bins] = hist(data, 100);
plot(bins, vals / trapz(bins, vals), '*', bins, 1/sqrt(2*pi)*exp(-bins.^2 /2), '-');
legend('normalInverseCDF', 'p(x)');
xlabel('x');
ylabel('density');
