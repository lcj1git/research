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

set title "RHF and FCI Energies of H_{10}, Two Rows"
set xlabel "{/Symbol a} (degrees)"
set ylabel "E (Hartree)"
set xtics 0, 5.0 
set grid xtics x2tics ytics
set key bottom left

set output "compareFCItoRHF.eps"
set xrange [45:85]
set yrange [-5.3:-4.6]
plot "h10_2rows_FCIsinglets" u 1:2 w l ls 1 ti "FCI n=1",  \
     "h10_2rows_FCIsinglets" u 1:3 w l ls 2 ti "FCI n=2",  \
     "h10_2rows_FCIsinglets" u 1:4 w l ls 3 ti "FCI n=3",  \
     "alpha25/h10RHFE" u 1:3 w l ls 4 ti "RHF 2A_g 3B_u",  \
     "alpha45/h10RHFE" u 1:3 w l ls 5 ti "RHF 3A_g 2B_u"
