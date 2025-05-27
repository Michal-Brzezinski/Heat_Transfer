#include "math.h"
#include "CALERF.h" 
#include "UTILS.h"


void utilspack::warunek_poczatkowy(long double* U, const long double* X, int N) {
    //-------------------------------------------------------------------
    // Warunek początkowy U(x,0):
    // Funkcja inicjalizuje wartości dla tablicy U na odpowiednie 
    // biorąc pod uwagę podany warunek początkowy
    //
    //  Argumenty:
    //      U       - tablica, w której zapisywane będą wartości początkowe
    //      X       - tablica przechowująca wartości węzłów przestrzennych
    //      N       - liczba węzłów siatki przestrzennej
    //
    //  Zwraca: Nic
    //-------------------------------------------------------------------

    for (int i = 0; i < N; i++) {
        U[i] = (X[i] < 0.0L) ? 0.0L : expl(-X[i] / b);
    }
}



long double utilspack::rozwiazanie_analityczne(long double x, long double t, int N) {
    //-------------------------------------------------------------------
    // Rozwiązanie analityczne równania dyfuzji:
    // znalezienie wartości U(x,t) dla podanych argumentów
    //
    //  Argumenty:
    //      x       - zmienna przestrzenna
    //      t       - zmienna czasu
    //      N       - liczba węzłów siatki przestrzennej
    //
    //  Zwraca: Wartość long double obliczonego analitycznie rozwiązania 
    //          dla podanego w treści zadania wzoru
    //-------------------------------------------------------------------

    long double z       = (2.0L * D * t / b - x) / (2.0L * sqrtl(D * t));

    long double pref    = 0.5L * expl(D * t / (b * b) - x / b);
    
    return pref * calerfpack::erfc_LD(z);
    //  używana jest funkcja z pakietu CALERF, udostępnionego przez prowadzącego
}



long double utilspack::compute_max_error(const long double* U_num, const long double* X, long double t, int N) {
    //-------------------------------------------------------------------
    //  Funkcja oblicza maksymalny błąd między rozwiązaniem 
    //  numerycznym a analitycznym w danym czasie t, obliczając
    //  wartosci bezwzgl. błędów dla każdego węzła siatki przestrzennej
    //  i zwracając największy z nich.
    //
    //  Argumenty:
    //      U_num   - Tablica przybliżonych wartości funkcji dla danego poziomu czasu
    //      X       - Tablica węzłów (siatka przestrzenna)
    //      t       - zadany poziom czasowy
    //      N       - liczba węzłów siatki przestrzennej
    //
    //  Zwraca:  MAKSYMALNY BŁĄD BEZWZGLĘDNY na danym poziomie czasowym
    //           pomiędzy obliczonym - przybliżonym rozwiązaniem na danym 
    //           poziomie czasowym i jego rozwiązaniem analutycznym.
    //-------------------------------------------------------------------

    long double max_err = 0.0L;
    for (int i = 0; i < N; ++i) {
        long double ue = rozwiazanie_analityczne(X[i], t, N);
        long double e = fabsl(U_num[i] - ue);
        if (e > max_err) {
            max_err = e;
        }
    }
    return max_err;
}