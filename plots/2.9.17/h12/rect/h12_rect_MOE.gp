set term postscript enhanced eps color "Palatino" 16
set size 0.7, 0.7
set encoding iso

# Set line styles
set style line  1 lt 1 lw 3 lc rgb "#0070FF" pt 6    # blue
set style line  2 lt 1 lw 3 lc rgb "#E52B50" pt 6    # red
set style line  3 lt 1 lw 3 lc rgb "#009A63" pt 6    # green
set style line  4 lt 1 lw 3 lc rgb "#702963" pt 6    # purple
set style line  5 lt 1 lw 3 lc rgb "#F28500" pt 6    # orange
set style line  6 lt 1 lw 3 lc rgb "#00BFFF" pt 6    # deep sky blue
set style line  7 lt 1 lw 3 lc rgb "#FF91A4" pt 6    # pink
set style line  8 lt 1 lw 3 lc rgb "#50C878" pt 6    # light green
set style line  9 lt 1 lw 3 lc rgb "#C9A0DC" pt 6    # lavender
set style line 10 lt 1 lw 3 lc rgb "#000000" pt 6    # black
set style line 11 lt 2 lw 3 lc rgb "#0070FF" pt 6    # blue
set style line 12 lt 2 lw 3 lc rgb "#E52B50" pt 6    # red
set style line 13 lt 2 lw 3 lc rgb "#009A63" pt 6    # green
set style line 14 lt 2 lw 3 lc rgb "#702963" pt 6    # purple
set style line 15 lt 2 lw 3 lc rgb "#F28500" pt 6    # orange
set style line 16 lt 2 lw 3 lc rgb "#00BFFF" pt 6    # deep sky blue
set style line 17 lt 2 lw 3 lc rgb "#FF91A4" pt 6    # pink
set style line 18 lt 2 lw 3 lc rgb "#50C878" pt 6    # light green
set style line 19 lt 2 lw 3 lc rgb "#C9A0DC" pt 6    # lavender
set style line 20 lt 2 lw 3 lc rgb "#000000" pt 6    # black
set style line 100 lt 1 lw 1 lc rgb "#FFFFFF"

set title "H_{12} Molecular Orbital Energies"
set xlabel "Alpha (degrees)"
set ylabel "E Hartree)"
set xtics 0, 10.0
#set x2tics ( 45.00 ) format "" scale 0
set grid xtics x2tics ytics

set output "h12_rect_MOE.eps"
set xrange [0:*]
set yrange [*:*]
plot "h12_rect_MOE" u 1:3 w l ls 1 ti "2^{nd}HOMO",  \
     "h12_rect_MOE" u 1:4 w l ls 2 ti "HOMO",  \
     "h12_rect_MOE" u 1:5 w l ls 3 ti "LUMO", \
     "h12_rect_MOE" u 1:6 w l ls 4 ti "2^{nd}LUMO"
