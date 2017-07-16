set xlabel 'N'
set ylabel 'absolute error'
set grid
set logscale x
set logscale y
set size ratio 1.0
set style data lines
plot "../output/testfunction_task15_error.txt" using 1:2 lt rgb "blue" title "CC RW" with linespoints, "../output/testfunction_task15_error.txt" using 1:3 lt rgb "green" title "CC BB" with linespoints
