#   Wywołanie: 
#   gnuplot "ML_Thomas_max_error_vs_t_plot.gp"



# Ustaw separator na przecinek (CSV)
set datafile separator ","

# Jeśli plik zawiera nagłówek, można ustawić:
set key autotitle columnhead

# Ustawienia osi i tytułu
set xlabel "t"
set ylabel "error_{max}"
set title "Thomas: Zaleznosc error_{max} od t"

# zakresy osi
#set xrange [-1:1]
#set xrange [-5:0]

set xtics -0, 0.2, 1.5
set ytics -1, 0.1, 1.5
set grid

# Ustawienie równej proporcji osi X i Y
set size ratio -1
set size square

set terminal qt size 600,600


# Rysujemy dane zapisane w pliku
plot "ML_Thomas_maxerror_vs_time.csv" using 1:2 with linespoints lw 2 pt 7 title "error(t)"

pause -1 "Nacisnij dowolny klawisz, aby zakonczyc"
