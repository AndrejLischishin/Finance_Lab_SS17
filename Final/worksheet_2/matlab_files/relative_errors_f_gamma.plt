set xlabel 'N'
set ylabel 'relative error'
set grid
set logscale x
set logscale y
set size ratio 1.0
set style data lines
plot "../output/relative_errors_f_gamma.txt" using 1:2 lt rgb "blue" title "Monte Carlo" with linespoints, "../output/relative_errors_f_gamma.txt" using 1:3 lt rgb "green" title "Trapezoidal rule" with linespoints, "../output/relative_errors_f_gamma.txt" using 1:4 lt rgb "red" title "Clenshaw Curtis" with linespoints, "../output/relative_errors_f_gamma.txt" using 1:5 lt rgb "yellow" title "Gauss Legendre" with linespoints
