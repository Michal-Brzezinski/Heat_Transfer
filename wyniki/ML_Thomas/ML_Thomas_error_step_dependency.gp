#   KOMENDA wywołania: 
#   gnuplot "ML_Thomas_error_step_dependency.gp"




# Ustaw separator na przecinek (CSV)
set datafile separator ","

# Jeśli plik zawiera nagłówek, można ustawić:
set key autotitle columnhead

# Ustawienia osi i tytułu
set xlabel "log_{10}(h)"
set ylabel "log_{10}(error_{max})"
set title "ML Thomas: Zaleznosc log_{10}(error_{max}) od log_{10}(h)"

# zakresy osi
#set xrange [-1:1]
#set xrange [-5:0]

set xtics -2, 0.5, 0
set ytics -5, 0.5, -1
set grid

# Ustawienie równej proporcji osi X i Y
set size ratio -1
#set size square

set terminal qt size 600,600


# Rysujemy dane zapisane w pliku
plot "ML_Thomas_results_error_step.csv" using 1:2 with linespoints lw 2 pt 7 title "Doswiadczalny rzad dokladnosci"

# Dopasowanie prostej: funkcja liniowa f(x)=A*x+B
f(x) = A*x + B
fit f(x) "ML_Thomas_results_error_step.csv" using 1:2 via A, B

# Nakładamy na wykres dopasowaną linię
replot f(x) title sprintf("Dopasowanie: y = %.3f x + %.3f", A, B)

pause -1 "Nacisnij dowolny klawisz, aby zakonczyc"
