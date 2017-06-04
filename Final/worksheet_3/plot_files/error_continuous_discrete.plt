set xlabel 'M'
set ylabel 'absolute error'
set grid
set logscale x
set logscale y
set size ratio 1.0
set style data lines
plot "../output/error_continuous_discrete.txt" using 1:2 lt rgb "blue" title "Error between discrete and continuous geometric average" with linespoints
