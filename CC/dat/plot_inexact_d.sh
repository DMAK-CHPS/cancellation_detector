set xrange [0:52]
set xlabel 'taille de la cancellation detectée'
set ylabel 'valeurs retournées par inexact (normalisées)' 
plot 'dat/inexact_d.dat' title 'inexact double autour de 1'
pause -1 "press enter to close"