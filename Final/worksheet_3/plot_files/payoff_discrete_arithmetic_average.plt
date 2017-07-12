set key off
set grid
set title "Payoff of discrete arithmetic average"
splot '../output/plot_payoff_discrete_arithmetic_average.txt' using 1:2:3 with points palette pointsize 0.5 pointtype 7
