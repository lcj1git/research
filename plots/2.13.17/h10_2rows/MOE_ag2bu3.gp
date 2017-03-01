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

set title "RHF MO Energies of H_{10}, Two Rows, a_g^2b_u^3"
set xlabel "{/Symbol a} (degrees)"
set ylabel "E (Hartree)"
set xtics 0, 5.0
set xtics add ("XC" 19.34) 
set xtics add ("XC" 37.44) 
set xtics add ("XC" 52.56) 
set xtics add ("XC" 70.66) 
set grid xtics x2tics ytics

set output "MOE_ag2bu3.eps"
set xrange [45:85]
set yrange [*:*]
plot "ag2bu3/h10MOE" u 1:6 w l ls 1 ti "2^{nd} HOMO",  \
     "ag2bu3/h10MOE" u 1:7 w l ls 2 ti "HOMO",  \
     "ag2bu3/h10MOE" u 1:8 w l ls 3 ti "LUMO",  \
     "ag2bu3/h10MOE" u 1:9 w l ls 4 ti "2^{nd} LUMO"