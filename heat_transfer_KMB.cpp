#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <utility> // dla std::swap

//  Pakiet udostęniony przez prowadzącego
#include "pakiety/CALERF.h" 

//  Pakiet dodatkowy (programu użytkowe)
#include "pakiety/UTILS.h"

/*  
    Komenda do kompilacji kodu: 
    g++ heat_transfer_KMB.cpp pakiety/calerf.cpp pakiety/utils.cpp -o KMB -lstdc++

    Komenda wykonująca program:
    ./KMB
*/

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
    
    // warunki brzegowe
    U_new[0] = 0.0L;
    U_new[N-1] = 0.0L;
    

    for (int i = 1; i + 1 < N; ++i) {
        U_new[i] = U_old[i] + lambda * (U_old[i + 1] - 2.0L * U_old[i] + U_old[i - 1]);
        //  Jest to przekształcony wzór KMB
    }

}

int main() {

    // Alokacja tablic dynamicznych
    long double* X   = new long double[N];  //  tablica przechowująca wartości węzłów siatki przestrzennej
    long double* U  = new long double[N];  //  tablica przechowująca wartości funkcji dla KMB
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
    //  Przykładem są wartości: M=1000, N=240

    // Wypisanie wymiarów siatki i lambdy
    std::cout << "N = " << N << ", M = " << M << ", lambda = " << lambda << std::endl;

    // Pętla czasowa
    for (int n = 0; n < M; ++n) {
        // Metoda KMB
        oblicz_nastepny_poziom_czasowy_KMB(U, Tmp, lambda);
        // Zamiana wskaźników, aby uniknąć kopiowania tablic – teraz Ue wskazuje na wynik nowej iteracji
        std::swap(U, Tmp);

    }

    // Obliczenie błędów przy czasie t_max
    long double err_kmb = utilspack::compute_max_error(U, X, t_max);
    std::cout << "Max error KMB = " << err_kmb << std::endl;

    // Zapis wyników do pliku CSV
    std::ofstream fout("KMBresults.csv");
    fout << "x,U_explicit,U_exact\n";
    for (int i = 0; i < N; ++i) {
        
        long double u_exact = utilspack::rozwiazanie_analityczne(X[i], t_max);
        fout << X[i] << "," << U[i] <<  "," << u_exact << "\n";
    }
    fout.close();

    // Dealokacja pamięci
    delete[] X;
    delete[] U;
    delete[] Tmp;

    return 0;
}
