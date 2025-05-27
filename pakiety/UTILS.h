#ifndef __utils_h
#define __utils_h

//----------------------------------------------------------------------
// Stałe fizyczne i przestrzenne
//----------------------------------------------------------------------
const long double D     = 1.0e0L;    // współczynnik dyfuzji
const long double b     = 1.0e-1L;   // parametr początkowy
const long double t_max = 1.0e0L;    // czas końcowy
const long double a     = 6.0e0L;   // a >= 6 * sqrt(D*t_max) lub większe


//----------------------------------------------------------------------
// Funkcje użytkowe we wszystkich podprogramach KMB i ML
//----------------------------------------------------------------------
namespace utilspack{
    void warunek_poczatkowy(long double* U, const long double* X, int N);
    long double compute_max_error(const long double* U_num, const long double* X, long double t, int N);
    long double rozwiazanie_analityczne(long double x, long double t, int N);
}

#endif