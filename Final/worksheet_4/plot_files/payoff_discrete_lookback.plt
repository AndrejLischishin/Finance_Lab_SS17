set key off
set grid
set title "Payoff of discrete Lookback option"
splot '../output/plot_integrand_lookback.txt' using 1:2:3 with points palette pointsize 0.5 pointtype 7
