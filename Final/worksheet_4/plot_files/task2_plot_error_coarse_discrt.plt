set xlabel 'N'
set ylabel 'error[|ref.val. - calc.res.|]'
set grid
set logscale x
set logscale y
set size ratio 1.0
set style data lines
set title "M = 64, barrier down-out call, QMC, MC, coarse discret."
plot "../output/task2_error_coarse_discretisation.txt" using 1:2 lt rgb "blue" title "QMC RW" with linespoints, "../output/task2_error_coarse_discretisation.txt" using 1:3 lt rgb "green" title "QMC BB" with linespoints, "../output/task2_error_coarse_discretisation.txt" using 1:4 lt rgb "red" title "MC RW" with linespoints, "../output/task2_error_coarse_discretisation.txt" using 1:5 lt rgb "black" title "QMC BB" with linespoints
