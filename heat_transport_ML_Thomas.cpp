#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <utility> // dla std::swap

//  Pakiet udostęniony przez prowadzącego
#include "pakiety/CALERF.h" 

//  Pakiet dodatkowy (programu użytkowe)
#include "pakiety/UTILS.h"

//  Pakiet dodatkowy (stworzony na zajęciach laboratoryjnych)
#include "pakiety/THOMAS.h"

// INFO !!! : DLA M=1000, N=380 jest absurdalnie duża różnica w porównaniu do LU dla pełnej macierzy

/*  
    Komenda do kompilacji kodu: 
    g++ heat_transport_ML_Thomas.cpp pakiety/calerf.cpp pakiety/utils.cpp pakiety/thomas.cpp -o ML_Thomas -lstdc++

    Komenda wykonująca program:
    ./ML_Thomas
*/

void oblicz_nastepny_poziom_czasowy_Laasonen_Thomas(const long double* U_old, long double* U_new, long double lambda) {
    //-------------------------------------------------------------------
    //  Funkcja oblicza przybliżoną wartość funkcji na kolejnym poziomie czasowym
    //  Oblicza układ równań z macierzą trójdiagonalną za pomocą:
    //      a) Algorytmu Thomasa
    //      b) Algorytmu dekompozycji LU dla macierzy pełnej
    //  Warunki brzegowe: U_new[0]=U_new[N-1]=0 (ustawiane na początku pętli for)
    
    //  Argumenty:
    //      U_old   - Tablica wartości funkcji dla bieżącego poziomu czasu
    //      U_new   - Tablica wartości funkcji dla nowego poziomu czasu
    //      lambda  - parametr lambda: D*dt/h^2

    //  Zwraca: Nic -> operacje na wskaźnikach
    //-------------------------------------------------------------------
    
    // Alokacja tablic na współczynniki układu trójdiagonalnego
    long double* l = new long double[N]; // dolna przekątna
    long double* d = new long double[N]; // główna przekątna
    long double* u = new long double[N]; // górna przekątna
    long double* c = new long double[N]; // wyrazy wolne
    

    for (int i = 0; i < N; ++i) {
        // Uzupełnienie macierzy A (a właściwie jej diagonali) odpowiednimi wyrazami

        if (i == 0 || i == N - 1) {
            //  Warunki brzegowe: U = 0 na brzegach(pierwszy i ostatni węzeł)
            //  Poniższe przekształcenie wynika bezpośrednio z postaci 
            //  macierzy A w metodzie Laasonen, gdzie 1. i ostatni wiersz
            //  odpowiadają za wartości funkcji na brzegach

            l[i] = 0.0L;
            d[i] = 1.0L;
            u[i] = 0.0L;
            c[i] = 0.0L;

        } else {
            //  Pozostałe wyrazy macierzy A są tutaj obliczane.
            //  Poniższe wynika z przekształcenia równania w metodzie Laasonen
            
            l[i] = -lambda;
            d[i] = 1.0L + 2.0L * lambda;
            u[i] = -lambda;
            c[i] = U_old[i];
        }
    }
    
    // Rozwiązujemy układ trójdiagonalny wykorzystując pakiet thomaspack:
    // Parametry: N, l, d, u, b (c), x (U_new)
    thomasldpack::Thomas(N, l, d, u, c, U_new);
    
    delete[] l;
    delete[] d;
    delete[] u;
    delete[] c;
    //  Zwolnienie zbędnych zasobów
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
        
        oblicz_nastepny_poziom_czasowy_Laasonen_Thomas(U, Tmp, lambda);
        std::swap(U, Tmp);
    }

    // Obliczenie błędów przy czasie t_max
    long double err_l = utilspack::compute_max_error(U, X, t_max);
    std::cout << "Max error laasonen Thomas = " << err_l << std::endl;

    // Zapis wyników do pliku CSV
    std::ofstream fout("ML_Thomas_results.csv");
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
