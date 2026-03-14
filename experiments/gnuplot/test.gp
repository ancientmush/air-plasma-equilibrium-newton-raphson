unset key
set xrange [0:20000]
set yrange [0:1]
set xlabel 'Temperature [K]'
set xlabel font "Times New Roman, 25"
set ylabel 'Mole Fraction [-]'
set ylabel font "Times New Roman, 25"
#set size square
set tics font "Times New Roman, 20"
set border linewidth 2
set ylabel offset -1,0
set bmargin 5
plot for [i=2:12] 'pressure_test.dat' using 1:i with lines lw 3 lt -1 linecolor rgb 'black'

