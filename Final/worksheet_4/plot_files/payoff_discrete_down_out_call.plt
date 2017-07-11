set key off
set grid
set title "Payoff of discrete Down-Out Call option"
splot '../output/plot_integrand_down_out_call.txt' using 1:2:3 with points palette pointsize 0.5 pointtype 7
