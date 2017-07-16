set xlabel 'N'
set ylabel 'error[|ref.val. - calc.res.|]'
set grid
set logscale x
set logscale y
set size ratio 1.0
set style data lines
set title "M = 64, lookback, QMC, MC"
plot "../output/task6_error.txt" using 1:2 lt rgb "blue" title "QMC BB" with linespoints, "../output/task6_error.txt" using 1:3 lt rgb "green" title "MC BB" with linespoints