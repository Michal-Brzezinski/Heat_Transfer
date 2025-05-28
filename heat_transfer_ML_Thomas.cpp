#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <utility> // dla std::swap
#include <set>
#include <string>

//  Pakiet udostęniony przez prowadzącego
#include "pakiety/CALERF.h" 

//  Pakiet dodatkowy (programu użytkowe)
#include "pakiety/UTILS.h"

//  Pakiet dodatkowy (stworzony na zajęciach laboratoryjnych)
#include "pakiety/THOMAS.h"


/*  
    Komenda do kompilacji kodu: 
    g++ heat_transfer_ML_Thomas.cpp "pakiety/CALERF.cpp" "pakiety/UTILS.cpp" "pakiety/THOMAS.cpp" -o ML_Thomas

    Komenda wykonująca program:
    ./ML_Thomas
*/


//___________________________________________________________________________________________________
//  WSTĘPNA KONFIGURACJA DLA PUNKTÓW 2 I 3
//  ODKOMENTOWAĆ DLA WYKONANIA PKT 2 I 3 :
#define POINT_2_AND_3

#ifdef POINT_2_AND_3
        //----------------------------------------------------------------------
        // Rozmiary siatki
        //----------------------------------------------------------------------

        //  liczba węzłów siatki przestrzennej
        int Xs = 2371;  //  np.: 240 -> KMB, 380 -> ML

        //  liczba węzłów siatki czasowej 
        int Ts = 39039; //  np. 1000 -> KMB, ML 

        //----------------------------------------------------------------------
        // Wartości kroków na siatce czasowo-przestrzennej
        //----------------------------------------------------------------------
        long double h  = static_cast<long double>((2.0L*a)/(Xs-1));     // krok przestrzenny
        long double dt = static_cast<long double>(t_max/(Ts-1));     // krok czasowy (krok całkowania)
#endif
//____________________________________________________________________________________________________



void oblicz_nastepny_poziom_czasowy_Laasonen_Thomas(const long double* U_old, long double* U_new, long double lambda, const int N) {
    //-------------------------------------------------------------------
    //  Funkcja oblicza przybliżoną wartość funkcji na kolejnym poziomie czasowym
    //  Oblicza układ równań z macierzą trójdiagonalną za pomocą:
    //      a) Algorytmu Thomasa
    //      b) Algorytmu dekompozycji LU dla macierzy pełnej
    //  Warunki brzegowe: U_new[0]=U_new[N-1]=0 (ustawiane na początku pętli for)
    //
    //  Argumenty:
    //      U_old   - Tablica wartości funkcji dla bieżącego poziomu czasu
    //      U_new   - Tablica wartości funkcji dla nowego poziomu czasu
    //      lambda  - parametr lambda: D*dt/h^2
    //      N - liczba węzłów siatki przestrzennej
    //
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
    thomaspack::Thomas(N, l, d, u, c, U_new);
    
    delete[] l;
    delete[] d;
    delete[] u;
    delete[] c;
    //  Zwolnienie zbędnych zasobów
}


#ifdef POINT_2_AND_3

int main() {

    // Alokacja tablic dynamicznych
    long double* T   = new long double[Ts];  //  tablica przechowująca wartości węzłów siatki czasowej
    long double* X   = new long double[Xs];  //  tablica przechowująca wartości węzłów siatki przestrzennej
    long double* U  = new long double[Xs];  //  tablica przechowująca wartości funkcji dla KMB
    long double* Tmp = new long double[Xs];  //  tablica przechowująca tymczasowe wartości funkcji
    int i = 0; //  zmienna iteracyjna, aby nie definiować ciągle nowej

    // Utworzenie siatki przestrzennej jako: X[i] = -a + i*h
    for (i = 0; i < Xs; ++i) {
        X[i] = -a + static_cast<long double>(i) * h;
    }

    // Utworzenie siatki czasowej jako: T[i] = i*dt
    for (i = 0; i < Ts; ++i) {
        T[i] = static_cast<long double>(i) * dt;
    }
    
    // Inicjalizacja warunku początkowego U(x,0)
    utilspack::warunek_poczatkowy(U, X, Xs);

    //  definiujemy parametr lambda dla obu metod
    long double lambda = D * dt / (h * h);
    //  Do obliczenia optymalnego lambda warto znać zależności: w ML (dla lambda bliskiego 1): 144M = N^2

    // Wypisanie wymiarów siatki i lambdy
    std::cout << "węzłów przestrzennych: " << Xs << ", węzłów czasowych: " << Ts << ", lambda = " << lambda << std::endl;


    std::set<int> save_indexes= {0, 1, 10, 30, 80, 1000, 5000, 10000, Ts-1};
    //  Tablica do przechowywania indeksów iteracji, w których zapisywane są wyniki

    
    // otwarcie pliku przed rozpoczęciem petli, aby nie nadpisywać pliku
    std::ofstream file_errr_time("wyniki/ML_Thomas/ML_Thomas_maxerror_vs_time.csv");
    file_errr_time << "t,e_max\n";
    long double err_kmb;
    // te instrukcje przed pętlą aby uniknąc redundancji danych w kodzie

    // Pętla czasowa
    for (int n = 0; n < Ts; n++) {
        // Metoda KMB
        std::string template_filename = "wyniki/ML_Thomas/ML_Thomas_results";
        if(save_indexes.count(n)){
            // Jeśli podany indeks jest jednym z wybranych do zapisu to zapisz do pliku CSV:

            //-------------------------- ZAPIS DO PLIKU CSV -------------------------------------
            std::ofstream fout(template_filename + std::to_string(n) + "iter.csv");   // np. LU_Thomas_results0iter.csv
            
            fout << "x,U_ML_Thomas,U_exact\n";
            for (i = 0; i < Xs; i++) {
                
                long double u_exact = utilspack::rozwiazanie_analityczne(X[i], T[n], Xs);
                fout << X[i] << "," << U[i] <<  "," << u_exact << "\n";
            }
            fout.close();
            //------------------------------------------------------------------------------------
        }
        
        //----------------- ZAPISANIE KROKU CAŁKOWANIA I BŁĘDU DO PLIKU CSV ------------------
        err_kmb = utilspack::compute_max_error(U, X, T[n], Xs);
        file_errr_time << T[n] << "," << err_kmb <<"\n";
        //------------------------------------------------------------------------------------

        oblicz_nastepny_poziom_czasowy_Laasonen_Thomas(U, Tmp, lambda, Xs);
        // Zamiana wskaźników, aby uniknąć kopiowania tablic – teraz Ue wskazuje na wynik nowej iteracji
        std::swap(U, Tmp);

    }
    file_errr_time.close();

    // Dealokacja pamięci
    delete[] T;
    delete[] X;
    delete[] U;
    delete[] Tmp;

    return 0;
}
#endif
