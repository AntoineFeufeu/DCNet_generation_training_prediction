 set term wxt
 #set term postscript landscape color solid "Helvetica" 22
 #set output "energy.ps"
 # set xrange [0:60]
 set logscale y
 set xlabel "Time (s)"
 set ylabel "Energy (J)"
 set loadpath "./OUTPUT_FILES/"
plot "energy.dat" us 1:4 t "Total Energy" w l lc 1, "energy.dat" us 1:3 t "Potential Energy" w l lc 2
 pause -1 "Hit any key..."
