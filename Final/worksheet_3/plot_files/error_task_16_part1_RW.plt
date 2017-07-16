set xlabel 'N'
set ylabel 'absolute error'
set grid
set logscale x
set logscale y
set size ratio 1.0
set style data lines
set title "All methods using RW, discrete geometric Asian option“
plot "../output/testfunction_task16_error.txt" using 1:2 lt rgb "blue" title "QMC" with linespoints, "../output/testfunction_task16_error.txt" using 1:4 lt rgb "green" title "MC" with linespoints, "../output/testfunction_task16_error.txt" using 1:6 lt rgb "red" title "FG TRAP" with linespoints, "../output/testfunction_task16_error.txt" using 1:8 lt rgb "yellow" title "FG CC" with linespoints, "../output/testfunction_task16_error.txt" using 1:10 lt rgb "cyan" title "SG TRAP" with linespoints, "../output/testfunction_task16_error.txt" using 1:12 lt rgb "black" title "SG CC" with linespoints
