set xrange [0:52]
set xlabel 'taille des cancellations generées'
set ylabel 'valeurs retournées par inexact' 
plot 'dat/double.dat' title 'detection de cancellation et bruitage par inexact de reel double precison'
pause -1 "press enter to close"