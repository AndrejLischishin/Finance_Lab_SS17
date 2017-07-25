n = 100.0 # anzahl bins
max = 6 # maximaler wert
min = -6 # minimaler wert
rows = 1000000.0 # anzahl messwerte
sigma = 1 # standardabweichung
mu = 0 # mittelrwer

# gaussverteilung:
gauss(x) = 1. / (sigma * sqrt(2 * pi)) * exp(-(x - mu)**2 / (2 * sigma**2))

# breite der bins
width = (max - min) / n

# histogrammfunktion für gnuplot
hist(x, width) = width * floor(x / width) + width / 2.0

# ein paar gnuplot-einstellungen
set xlabel "x"
set ylabel "density"

plot '../output/rejection_sampl.txt' using (hist($1, width)):(100.0 / rows) smooth freq lc rgb"red" notitle, \
gauss(x) * width * 100.0 notitle lc rgb"blue"
