set xlabel 'N'
set ylabel 'error[|ref.val. - calc.res.|]'
set grid
set logscale x
set logscale y
set size ratio 1.0
set style data lines
set title "M=64, asian call, arithmetic, control variates, QMC, MC"
plot "../output/task7_error.txt" using 1:2 lt rgb "blue" title "QMC BB" with linespoints, "../output/task7_error.txt" using 1:3 lt rgb "green" title "QMC BB CV" with linespoints, "../output/task7_error.txt" using 1:4 lt rgb "red" title "MC BB" with linespoints, "../output/task7_error.txt" using 1:5 lt rgb "black" title "QMC BB CV" with linespoints
