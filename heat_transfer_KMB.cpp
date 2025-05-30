#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <utility> // dla std::swap
#include <set>  // dla zapisywania w csv (sprawdzania indeksów)
#include <string>
#include <iomanip>
//  Pakiet udostęniony przez prowadzącego
#include "pakiety/CALERF.h" 
//  Pakiet dodatkowy (programu użytkowe)
#include "pakiety/UTILS.h"

/*  
            Komenda do kompilacji kodu: 
            g++ heat_transfer_KMB.cpp pakiety/CALERF.cpp pakiety/UTILS.cpp -o KMB

            Komenda wykonująca program:
            ./KMB
*/


//______________________________________________________
//------------------------------------------------------
//#define POINT_1 // ZAKOMENTOWAĆ DLA WYKONANIA PKT 2 I 3 
//------------------------------------------------------
//______________________________________________________


//___________________________________________________________________________________________________
//  WSTĘPNA KONFIGURACJA DLA PUNKTÓW 2 I 3
//  
#ifndef POINT_1
        //----------------------------------------------------------------------
        // Rozmiary siatki
        //----------------------------------------------------------------------

        //  liczba węzłów siatki przestrzennej
        int Xs = 1500;  //  np.: 240 -> KMB, 380 -> ML

        //  liczba węzłów siatki czasowej 
        int Ts = 39063; //  np. 1000 -> KMB, ML 

        //----------------------------------------------------------------------
        // Wartości kroków na siatce czasowo-przestrzennej
        //----------------------------------------------------------------------
        long double h  = static_cast<long double>((2.0L*a)/(Xs-1));     // krok przestrzenny
        long double dt = static_cast<long double>(t_max/(Ts-1));     // krok czasowy (krok całkowania)
#endif
//____________________________________________________________________________________________________





void oblicz_nastepny_poziom_czasowy_KMB(const long double* U_old, long double* U_new, long double lambda, const int N) {
    //-------------------------------------------------------------------
    // Funkcja oblicza przybliżoną wartość funkcji na kolejnym poziomie czasowym
    // Warunki brzegowe: U_new[0]=U_new[N-1]=0 (przyjmujemy, że już są ustawione)
    
    //  Argumenty:
    //  U_old - Tablica wartości funkcji dla bieżącego poziomu czasu
    //  U_new - Tablica wartości funkcji dla nowego poziomu czasu
    //  lambda - parametr lambda: D*dt/h^2
    //  N - liczba węzłów siatki przestrzennej

    //  Zwraca: Nic -> operacje na wskaźnikach
    //-------------------------------------------------------------------
    
    // warunki brzegowe
    U_new[0] = 0.0L;
    U_new[N-1] = 0.0L;
    

    for (int i = 1; i + 1 < N; ++i) {
        U_new[i] = U_old[i] + lambda * (U_old[i + 1] - 2.0L * U_old[i] + U_old[i - 1]);
        //  Jest to przekształcony wzór KMB
    }

}


#ifdef POINT_1

int main() {
    int Xs, Ts;  // zmienne przechowujące ilość węzłów
    long double h, dt;  // kroki

    std::ofstream fout("wyniki/KMB/KMBresults_error_step.csv");
    fout << "log10(h),log10(max_error)\n";
    fout << std::fixed << std::setprecision(19);

    for (int k = 1; k <= 50; ++k) {
        Xs = 24 * k;          // N jako wielokrotność 24
        Ts = 10 * k * k;       // M z zależności 576M = 10N^2

        // Używamy std::vector do alokacji pamięci
        
        long double* X   = new long double[Xs];  //  tablica przechowująca wartości węzłów siatki przestrzennej
        long double* U  = new long double[Xs];  //  tablica przechowująca wartości funkcji dla KMB
        long double* Tmp = new long double[Xs];  //  tablica przechowująca tymczasowe wartości funkcji

        h = (2.0L * a) / (Xs - 1);      // krok przestrzenny
        dt = t_max / (Ts - 1);          // krok czasowy (krok całkowania)

        // Utworzenie siatki przestrzennej: X[i] = -a + i*h
        for (int i = 0; i < Xs; ++i) {
            X[i] = -a + static_cast<long double>(i) * h;
        }

        // Inicjalizacja warunku początkowego U(x,0)
        utilspack::warunek_poczatkowy(U, X, Xs);

        // Parametr lambda dla metody
        long double lambda = D * dt / (h * h);

        // Wypisanie wymiarów siatki i lambdy
        std::cout << "węzłów przestrzennych: " << Xs << ", węzłów czasowych: " << Ts << ", lambda = " << lambda << std::endl;

        // Pętla czasowa (Ms-1 kroków)
        for (int n = 0; n < Ts - 1; ++n) {
            oblicz_nastepny_poziom_czasowy_KMB(U, Tmp, lambda, Xs);
            std::swap(U, Tmp);
}

        // Obliczenie błędu przy czasie t_max
        long double err_kmb = utilspack::compute_max_error(U, X, t_max, Xs);
        std::cout << "Max error KMB = " << err_kmb << std::endl;
        fout << log10(h) << "," << log10(err_kmb) << "\n";

        delete[] X;
        delete[] U;
        delete[] Tmp;

    }

    fout.close();
    
    return 0;
}
#endif



#ifndef POINT_1

int main() {

    auto start = std::chrono::high_resolution_clock::now(); // Zapisujemy czas rozpoczęcia

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
    //  Do obliczenia optymalnego lambda warto znać zależności: 
    //  (w KMB bliskiego 0.4) 576Ts = 10Xs^2 
    //  Przykładem są wartości: Ts=1000, Xs=240

    // Wypisanie wymiarów siatki i lambdy
    std::cout << "węzłów przestrzennych: " << Xs << ", węzłów czasowych: " << Ts << ", lambda = " << lambda << std::endl;


    std::set<int> save_indexes= {0, 1, 10, 30, 80, 200, 1000, 10000, Ts-1};
    //  Tablica do przechowywania indeksów iteracji, w których zapisywane są wyniki

    
    // otwarcie pliku przed rozpoczęciem petli, aby nie nadpisywać pliku
    std::ofstream file_errr_time("wyniki/KMB/KMB_maxerror_vs_time.csv");
    file_errr_time << "t,e_max\n";
    long double err_kmb;
    // te instrukcje przed pętlą aby uniknąc redundancji danych w kodzie

    // Pętla czasowa
    for (int n = 0; n < Ts; n++) {
        // Metoda KMB
        std::string template_filename = "wyniki/KMB/KMBresults";
        if(save_indexes.count(n)){
            // Jeśli podany indeks jest jednym z wybranych do zapisu to zapisz do pliku CSV:

            //-------------------------- ZAPIS DO PLIKU CSV -------------------------------------
            std::ofstream fout(template_filename + std::to_string(n) + "iter.csv");   // np. KMBresults0.csv
            
            fout << "x,U_KMB,U_exact\n";
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

        oblicz_nastepny_poziom_czasowy_KMB(U, Tmp, lambda, Xs);
        // Zamiana wskaźników, aby uniknąć kopiowania tablic – teraz Ue wskazuje na wynik nowej iteracji
        std::swap(U, Tmp);

    }
    file_errr_time.close();

    // Dealokacja pamięci
    delete[] T;
    delete[] X;
    delete[] U;
    delete[] Tmp;

    auto end = std::chrono::high_resolution_clock::now(); // Zapisujemy czas zakończenia
    std::chrono::duration<double> duration = end - start; // Obliczamy różnicę czasu
    std::cout << "Czas wykonania: " << duration.count() << " sekund\n";

    return 0;
}
#endif
