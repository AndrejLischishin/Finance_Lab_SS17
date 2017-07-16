set xlabel 'N'
set ylabel 'absolute error'
set grid
set logscale x
set logscale y
set size ratio 1.0
set style data lines
set title "discrete geometric Asian option, M=64"
plot "../output/testfunction_task17_error.txt" using 1:2 lt rgb "blue" title "QMC RW" with linespoints, "../output/testfunction_task17_error.txt" using 1:3 lt rgb "green" title "QMC BB" with linespoints, "../output/testfunction_task17_error.txt" using 1:4 lt rgb "red" title "MC RW" with linespoints, "../output/testfunction_task17_error.txt" using 1:5 lt rgb "black" title "MC BB" with linespoints
