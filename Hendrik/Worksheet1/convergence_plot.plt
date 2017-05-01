set xlabel 'N'
set ylabel 'estimation error'
set grid
set logscale x
set logscale y
set size ratio 1.0
set style data lines
plot "error1.txt" title "sigma = 0.1", "error2.txt" title "sigma = 2.0", "error3.txt" title "sigma = 10.0"
