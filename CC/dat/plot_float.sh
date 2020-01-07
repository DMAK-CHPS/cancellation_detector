set xrange [0:23]
set xlabel 'taille des cancellations generées'
set ylabel 'valeurs retournées par inexact' 
plot 'dat/float.dat' title 'detection de cancellation et bruitage par inexact de reel simple precison'
pause -1 "press enter to close"