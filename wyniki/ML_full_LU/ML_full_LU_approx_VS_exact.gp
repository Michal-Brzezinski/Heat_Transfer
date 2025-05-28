#   Wywołanie: 
#   gnuplot "ML_full_LU_approx_VS_exact.gp"



# Ustawienie separatora pól na przecinek
set datafile separator ","

# Jeśli plik zawiera nagłówek, można ustawić:
set key autotitle columnhead

# Ustawienia etykiet osi oraz legendy
set xlabel "x"
set ylabel "U"
set title "Porownanie przyblizen U_i z dokladnymi wartosciami U(x_i)"

set yrange [-1:6]
set grid

# Ustawienie równej proporcji osi X i Y
set size ratio -1
#set size square

set terminal qt size 600,600

# Wykres: U_KMB jako punkty, U_exact jako linia; używamy różnych kolorów
plot "ML_full_LU_results999iter.csv" using 1:2 with points pointtype 7 linecolor rgb "red" title "U_{KMB}", \
     "ML_full_LU_results999iter.csv" using 1:3 with lines linewidth 2 linecolor rgb "blue" title "U_{exact}"

pause -1 "Nacisnij dowolny klawisz, aby zakonczyc"