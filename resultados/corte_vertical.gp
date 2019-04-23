# Usar:
#  gnuplot -persists < corte_vertical.gp

#set term png               # (*)  Comenta para ver. Descomenta para guardar !!
#set output 'archivo.png'   # (**) Comenta para ver. Descomenta para guardar !!

set xlabel 'y'
set ylabel 'u1'
set xrange [0:1]
set yrange [-0.5:1]
set key top left

plot 'corte.txt'  using 1:2 title 'FEM Re=400' with lines ,\
     'ghia_ghia_shin_corte_vertical_Y-vs-U1.dat' using 1:3 title 'Ghia Re=400'
