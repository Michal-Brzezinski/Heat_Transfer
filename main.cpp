#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <utility> // dla std::swap

//  Pakiet udostęniony przez prowadzącego
#include "CALERF.h" 
//  Pakiet dodatkowy (stworzony na zajęciach laboratoryjnych)
#include "THOMAS.h"


/*  
    Komenda do kompilacji kodu: 
    g++ main.cpp calerf.cpp -o main -lstdc++  

    Komenda wykonująca program:
    ./main
*/


//----------------------------------------------------------------------
// Stałe fizyczne i przestrzenne
//----------------------------------------------------------------------
const long double D     = 1.0e0L;    // współczynnik dyfuzji
const long double b     = 1.0e-1L;   // parametr początkowy
const long double t_max = 1.0e0L;    // czas końcowy
// a = 6 * sqrt(D*t_max)
const long double a     = 6.0e0L;

//----------------------------------------------------------------------
// Rozmiary siatki
//----------------------------------------------------------------------
//  liczba węzłów siatki przestrzennej
const int N = 1200;

//  liczba węzłów siatki czasowej 
const int M = 1000;

//----------------------------------------------------------------------
// Wartości kroków na siatce czasowo-przestrzennej
//----------------------------------------------------------------------
const long double h  = (2*a)/(N-1);     // krok przestrzenny
const long double dt = t_max/(M-1);     // krok czasowy (krok całkowania)



//----------------------------------------------------------------------
// Funkcje pomocnicze
//----------------------------------------------------------------------


void warunek_poczatkowy(long double* U, const long double* X) {
    //-------------------------------------------------------------------
    // Warunek początkowy U(x,0):
    // Funkcja inicjalizuje wartości dla tablicy U na odpowiednie 
    // biorąc pod uwagę podany warunek początkowy

    //  Argumenty:
    //  U - tablica, w której zapisywane będą wartości początkowe
    //  X - tablica przechowująca wartości węzłów przestrzennych

    //  Zwraca: Nic
    //-------------------------------------------------------------------

    for (int i = 0; i < N; ++i) {
        U[i] = (X[i] < 0.0L) ? 0.0L : std::expl(-X[i] / b);
    }
}


long double rozwiazanie_analityczne(long double x, long double t) {
    //-------------------------------------------------------------------
    // Rozwiązanie analityczne równania dyfuzji:
    // znalezienie wartości U(x,t) dla podanych argumentów
    
    //  Argumenty:
    //  x - zmienna przestrzenna
    //  t - zmienna czasu

    //  Zwraca: Wartość long double obliczonego analitycznie rozwiązania 
    //          dla podanego w treści zadania wzoru
    //-------------------------------------------------------------------

    long double z       = (2.0L * D * t / b - x) / (2.0L * std::sqrtl(D * t));

    long double pref    = 0.5L * std::expl(D * t / (b * b) - x / b);
    
    return pref * calerfpack::erfc_LD(z);
    //  używana jest funkcja z pakietu CALERF, udostępnionego przez prowadzącego
}



void oblicz_nastepny_poziom_czasowy_KMB(const long double* U_old, long double* U_new, long double lambda) {
    //-------------------------------------------------------------------
    // Funkcja oblicza przybliżoną wartość funkcji na kolejnym poziomie czasowym
    // Warunki brzegowe: U_new[0]=U_new[N-1]=0 (przyjmujemy, że już są ustawione)
    
    //  Argumenty:
    //  U_old - Tablica wartości funkcji dla bieżącego poziomu czasu
    //  U_new - Tablica wartości funkcji dla nowego poziomu czasu
    //  lambda - parametr lambda: D*dt/h^2

    //  Zwraca: Nic -> operacje na wskaźnikach
    //-------------------------------------------------------------------

    for (int i = 1; i + 1 < N; ++i) {
        U_new[i] = U_old[i] + lambda * (U_old[i + 1] - 2.0L * U_old[i] + U_old[i - 1]);
        //  Jest to przekształcony wzór KMB
    }
}


void oblicz_nastepny_poziom_czasowy_Laasonen(const long double* U_old, long double* U_new, long double lambda) {
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

        if (i == 0 || i + 1 == N) {
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
    
    // Rozwiązujemy układ trójdiagonalny wykorzystując pakiet thomasldpack:
    // Parametry: N, l, d, u, b (c), x (U_new)
    thomaspack::Thomas(N, l, d, u, c, U_new);
    
    delete[] l;
    delete[] d;
    delete[] u;
    delete[] c;
    //  Zwolnienie zbędnych zasobów
}


long double compute_max_error(const long double* U_num, const long double* X, long double t) {
    //-------------------------------------------------------------------
    //  Funkcja oblicza maksymalny błąd między rozwiązaniem numerycznym a analitycznym w danym czasie t
    //
    //  Argumenty:
    //      U_num   - Tablica przybliżonych wartości funkcji dla danego poziomu czasu
    //      X       - Tablica węzłów (siatka przestrzenna)
    //      t       - zadany poziom czasowy
    //
    //  Zwraca:  MAKSYMALNY BŁĄD BEZWZGLĘDNY na danym poziomie czasowym
    //           pomiędzy obliczonym - przybliżonym rozwiązaniem na danym 
    //           poziomie czasowym i jego rozwiązaniem analutycznym.
    //-------------------------------------------------------------------

    long double max_err = 0.0L;
    for (int i = 0; i < N; ++i) {
        long double ue = rozwiazanie_analityczne(X[i], t);
        long double e = std::fabsl(U_num[i] - ue);
        if (e > max_err) {
            max_err = e;
        }
    }
    return max_err;
}



int main() {
    // Wypisanie wymiarów siatki
    std::cout << "N = " << N << ", M = " << M << std::endl;

    // Alokacja tablic dynamicznych
    long double* X   = new long double[N];  //  tablica przechowująca wartości węzłów siatki przestrzennej
    long double* U0  = new long double[N];  //  tablica przechowująca wartości funkcji dla warunku początkowego
    long double* Ue  = new long double[N];  //  tablica przechowująca wartości funkcji dla KMB
    long double* Ul  = new long double[N];  //  tablica przechowująca wartości funkcji  dla pośr. metody Laasonen
    long double* Tmp = new long double[N];  //  tablica przechowująca tymczasowe wartości funkcji

    // Utworzenie siatki przestrzennej jako: X[i] = -a + i*h
    for (int i = 0; i < N; ++i) {
        X[i] = -a + static_cast<long double>(i) * h;
    }
    
    // Inicjalizacja warunku początkowego U(x,0)
    warunek_poczatkowy(U0, X);

    // Kopiujemy U0 do Ue oraz Ul
    for (int i = 0; i < N; ++i) {
        Ue[i] = U0[i];
        Ul[i] = U0[i];
    }


    //  definiujemy parametr lambda dla obu metod
    long double lambda = D * dt / (h * h);
    //  Do obliczenia optymalnego lambda warto znać zależności: 
    //  (w KMB bliskiego 0.4) 576M = 10N^2
    //  Natomiast  w ML (dla lambda bliskiego 1): 144M = N^2

    // Pętla czasowa
    for (int n = 0; n < M; ++n) {
        // Metoda KMB
        oblicz_nastepny_poziom_czasowy_KMB(Ue, Tmp, lambda);
        // Zamiana wskaźników, aby uniknąć kopiowania tablic – teraz Ue wskazuje na wynik nowej iteracji
        std::swap(Ue, Tmp);
        
        // Metoda pośrednia Laasonen
        oblicz_nastepny_poziom_czasowy_Laasonen(Ul, Tmp, lambda);
        std::swap(Ul, Tmp);
    }

    // Obliczenie błędów przy czasie t_max
    long double err_e = compute_max_error(Ue, X, t_max);
    long double err_l = compute_max_error(Ul, X, t_max);
    std::cout << "Max error KMB = " << err_e << std::endl;
    std::cout << "Max error laasonen = " << err_l << std::endl;

    // Zapis wyników do pliku CSV
    std::ofstream fout("results.csv");
    fout << "x,U_explicit,U_laasonen,U_exact\n";
    for (int i = 0; i < N; ++i) {
        long double ue_x = rozwiazanie_analityczne(X[i], t_max);
        fout << X[i] << "," << Ue[i] << "," << Ul[i] << "," << ue_x << "\n";
    }
    fout.close();

    // Dealokacja pamięci
    delete[] X;
    delete[] U0;
    delete[] Ue;
    delete[] Ul;
    delete[] Tmp;

    return 0;
}
