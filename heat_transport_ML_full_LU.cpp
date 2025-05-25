#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <utility> // dla std::swap
#include "string.h" // do memset

//  Pakiet udostęniony przez prowadzącego
#include "pakiety/CALERF.h" 

//  Pakiet dodatkowy (programu użytkowe)
#include "pakiety/UTILS.h"

//  Pakiet dodatkowy (stworzony na zajęciach laboratoryjnych)
#include "pakiety/LU.h"

// INFO !!! : DLA M=1000, N=380, MIELI SIĘ Z 5 MINUT

/*  
    Komenda do kompilacji kodu: 
    g++ heat_transport_ML_full_LU.cpp pakiety/calerf.cpp pakiety/utils.cpp pakiety/lu.cpp -o ML_LU -lstdc++

    Komenda wykonująca program:
    ./ML_LU
*/

void oblicz_nastepny_poziom_czasowy_Laasonen_LU(const long double* U_old, 
                                              long double* U_new, 
                                              long double lambda) {
    //-------------------------------------------------------------------
    // Funkcja oblicza przybliżoną wartość funkcji na kolejnym poziomie czasowym
    // Metoda Laasonen – układ równań z macierzą trójdiagonalną (przy brzegach
    // ustawiamy U=0). W tej wersji budujemy pełną macierz A (w formacie jednowymiarowym)
    // i wektor prawej strony U_new (nasze per se "b"):
    //
    //   Dla wierszy brzegowych (i == 0 lub i == N-1):
    //       A[i, i] = 1, pozostałe elementy = 0, U_new[i] = 0.
    //
    //   Dla wierszy wewnętrznych (1 <= i <= N-2):
    //       A[i, i-1] = -lambda, A[i, i] = 1 + 2*lambda, A[i, i+1] = -lambda,
    //       U_new[i] = U_old[i]
    //
    // Po rozwiązaniu A*x = b, wynik (x) zostaje zapisany do U_new.
    //
    // Argumenty:
    //   U_old - wektor wartości funkcji dla bieżącego poziomu czasowego,
    //   U_new - wektor, do którego zapiszemy wynik kolejnej iteracji,
    //   lambda - parametr lambda: D*dt/h^2 (najlepiej bliski 1 dla tej metody)
    //-------------------------------------------------------------------
    
    // Alokujemy pełną macierz A o rozmiarze N x N (zapisywaną jako jednowymiarowa tablica):
    long double* A = new long double[N * N];
    
    // Zerujemy macierz A:
    memset(A, 0, N * N * sizeof(long double));

    
    // Budujemy macierz A oraz wektor b_vec:
    for (int i = 0; i < N; ++i) {
        
        if (i == 0 || i == (N - 1)) {
            //  Warunki brzegowe: U = 0 na brzegach(pierwszy i ostatni węzeł)
            //  Poniższe przekształcenie wynika bezpośrednio z postaci 
            //  macierzy A w metodzie Laasonen, gdzie 1. i ostatni wiersz
            //  odpowiadają za wartości funkcji na brzegach
            
            A[i * N + i] = 1.0L;
            U_new[i] = 0.0L;
        
        
        } else {
            // Wiersz i (wewnętrzny):
            A[i * N + (i - 1)] = -lambda;
            A[i * N + i]       = 1.0L + 2.0L * lambda;
            A[i * N + (i + 1)] = -lambda;
            U_new[i] = U_old[i];
        }
    }
    
    // Rozwiązujemy układ A*x = b_vec wykorzystując funkcję LU_decompose_and_solve.
    // Funkcja ta przyjmuje macierz A, wektor U_new (jako b) i rozmiar macierzy.
    ludecomposepack::LU_decompose_and_solve(A, U_new, N);

    // Zwolnienie alokowanych zasobów:
    delete[] A;
}


int main() {

    // Alokacja tablic dynamicznych
    long double* X   = new long double[N];  //  tablica przechowująca wartości węzłów siatki przestrzennej
    long double* U  = new long double[N];  //  tablica przechowująca wartości funkcji  dla pośr. metody Laasonen
    long double* Tmp = new long double[N];  //  tablica przechowująca tymczasowe wartości funkcji

    // Utworzenie siatki przestrzennej jako: X[i] = -a + i*h
    for (int i = 0; i < N; ++i) {
        X[i] = -a + static_cast<long double>(i) * h;
    }
    
    // Inicjalizacja warunku początkowego U(x,0)
    utilspack::warunek_poczatkowy(U, X);


    //  definiujemy parametr lambda dla obu metod
    long double lambda = D * dt / (h * h);
    //  Do obliczenia optymalnego lambda warto znać zależności: 
    //  (w KMB bliskiego 0.4) 576M = 10N^2 
    //  Natomiast  w ML (dla lambda bliskiego 1): 144M = N^2

    // Wypisanie wymiarów siatki i lambdy
    std::cout << "N = " << N << ", M = " << M << ", lambda = " << lambda << std::endl;

    // Pętla czasowa
    for (int n = 0; n < M; ++n) {
        // Metoda pośrednia Laasonen - obliczanie w pętli kolejnych wartości przybliżonych
        
        oblicz_nastepny_poziom_czasowy_Laasonen_LU(U, Tmp, lambda);
        std::swap(U, Tmp);
    }

    // Obliczenie błędów przy czasie t_max
    long double err_l = utilspack::compute_max_error(U, X, t_max);
    std::cout << "Max error laasonen LU = " << err_l << std::endl;

    // Zapis wyników do pliku CSV
    std::ofstream fout("ML_LUresults.csv");
    fout << "x,U_laasonen,U_exact\n";
    
    for (int i = 0; i < N; ++i) {
    
        long double u_exact = utilspack::rozwiazanie_analityczne(X[i], t_max);
        fout << X[i] << "," << U[i] << "," << u_exact << "\n";
    }
    fout.close();

    // Dealokacja pamięci
    delete[] X;
    delete[] U;
    delete[] Tmp;

    return 0;
}
