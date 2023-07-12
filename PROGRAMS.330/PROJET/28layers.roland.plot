set term X11
#set term postscript 
#set output "28layers.tdspdep.ps"
set xlabel "frequency (Hz)"
set ylabel "velocity (m/s)"
set yrange [0. : 4000.]
set xrange [0. : 250.]
plot '28layers.rc.dat' us ($2):(1000*$3) w l #with lines
