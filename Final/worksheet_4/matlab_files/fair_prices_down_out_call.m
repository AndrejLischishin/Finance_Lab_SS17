fid = fopen('../output/fair_prices_down_out_call.txt', 'r');
format_spec = '%f %f';
sizeA = [2 Inf];
A = fscanf(fid, format_spec, sizeA);
fclose(fid);
d = A(1,:);
x = A(2,:);

plot(d, x);

title('Fair prices for Down-Out call option');
ylabel('Price');
xlabel('B');
