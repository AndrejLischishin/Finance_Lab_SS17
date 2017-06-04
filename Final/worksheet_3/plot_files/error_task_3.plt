set xlabel 'N'
set ylabel 'absolute error'
set grid
set logscale x
set logscale y
set size ratio 1.0
set style data lines
plot "../output/error_task_3_M_10.txt" using 1:2 lt rgb "blue" title "M=10" with linespoints, "../output/error_task_3_M_200.txt" using 1:2 lt rgb "green" title "M=200" with linespoints
