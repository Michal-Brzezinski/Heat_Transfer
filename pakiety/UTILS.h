#ifndef __utils_h
#define __utils_h

//----------------------------------------------------------------------
// Stałe fizyczne i przestrzenne
//----------------------------------------------------------------------
const long double D     = 1.0e0L;    // współczynnik dyfuzji
const long double b     = 1.0e-1L;   // parametr początkowy
const long double t_max = 1.0e0L;    // czas końcowy
const long double a     = 6.0e0L;   // a = 6 * sqrt(D*t_max) lub większe

//----------------------------------------------------------------------
// Rozmiary siatki
//----------------------------------------------------------------------
//  liczba węzłów siatki przestrzennej
const int N = 380;  //  240 -> KMB, 380 -> ML

//  liczba węzłów siatki czasowej 
const int M = 1000; //  1000 -> KMB, ML 

//----------------------------------------------------------------------
// Wartości kroków na siatce czasowo-przestrzennej
//----------------------------------------------------------------------
const long double h  = static_cast<long double>((2.0L*a)/(N-1));     // krok przestrzenny
const long double dt = static_cast<long double>(t_max/(M-1));     // krok czasowy (krok całkowania)


//----------------------------------------------------------------------
// Funkcje użytkowe we wszystkich podprogramach KMB i ML
//----------------------------------------------------------------------
namespace utilspack{
    void warunek_poczatkowy(long double* U, const long double* X);
    long double compute_max_error(const long double* U_num, const long double* X, long double t);
    long double rozwiazanie_analityczne(long double x, long double t);
}

#endif