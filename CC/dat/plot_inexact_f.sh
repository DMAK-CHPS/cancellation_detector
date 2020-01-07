set xrange [0:23]
set xlabel 'taille de la cancellation detectée'
set ylabel 'valeurs retournées par inexact (normalisées)'
plot 'dat/inexact_f.dat' title 'inexact float autour de 1'
pause -1 "press enter to close"