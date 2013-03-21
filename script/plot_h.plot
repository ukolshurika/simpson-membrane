set term postscript eps color blacktext "Helvetica" 24
set output 'img/h.eps'
set xlabel 'h'
set ylabel 't'
set xr [1000000: 10000000000.0]
plot './data/h_last.out' with lines
set output
quit