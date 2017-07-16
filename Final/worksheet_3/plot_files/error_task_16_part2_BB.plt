set xlabel 'N'
set ylabel 'absolute error'
set grid
set logscale x
set logscale y
set size ratio 1.0
set style data lines
set title "All methods using BB, discrete geometric Asian option"
plot "../output/testfunction_task16_error.txt" using 1:3 lt rgb "blue" title "QMC" with linespoints, "../output/testfunction_task16_error.txt" using 1:5 lt rgb "green" title "MC" with linespoints, "../output/testfunction_task16_error.txt" using 1:7 lt rgb "red" title "FG TRAP" with linespoints, "../output/testfunction_task16_error.txt" using 1:9 lt rgb "yellow" title "FG CC" with linespoints, "../output/testfunction_task16_error.txt" using 1:11 lt rgb "cyan" title "SG TRAP" with linespoints, "../output/testfunction_task16_error.txt" using 1:13 lt rgb "black" title "SG CC" with linespoints
