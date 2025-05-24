#include <iostream>
#include <iomanip>
#include "ThomasLD.h"

using namespace std;



void thomasldpack::thomas_procedure_1(int N, const long double l[], long double d[], const long double u[]) {
    //-------------------------------------------------------------------
    //  Procedura operująca na macierzy A (trójdiagonalnej).
    //  W rzeczwywistości operacyjnej mamy 3 oddzielne tablice l,d,u
    //  Procedura dokonuje eliminacji w przód modyfikując główną przekątną.
    //  Dla i = 1,...,N-1: d[i] = d[i] - (l[i] / d[i-1]) * u[i-1]

    //  Argumenty:
    //  N - rozmiar macierzy A
    //  l[] - tablica wartości dolnej przekątnej macierzy 
    //  d[] - tablica wartości głównej przekątnej macierzy 
    //  u[] - tablica wartości górnej przekątnej macierzy 

    //  Zwraca: Nic -> operuje na wskaźnikach
    //-------------------------------------------------------------------
    
    for (int i = 1; i < N; i++) {
        long double m = l[i] / d[i - 1];
        d[i] = d[i] - m * u[i - 1];
    }
}



void thomasldpack::thomas_procedure_2(int N, const long double l[], const long double u[], 
        const long double d[], long double b[], long double x[]) {
    //-------------------------------------------------------------------
    // Procedura operująca na wektorze b.
    // Najpierw wykonuje eliminację w przód:
    // dla i = 1,...,N-1: b[i] = b[i] - (l[i] / d[i-1]) * b[i-1]
    // Następnie przeprowadza etap podstawiania wstecznego:
    // x[N-1] = b[N-1] / d[N-1];
    // dla i = N-2,...,0: x[i] = (b[i] - u[i]*x[i+1]) / d[i]

    //  Argumenty:
    //  N - rozmiar macierzy A -> zatem także wektora b
    //  l[] - tablica wartości dolnej przekątnej macierzy 
    //  u[] - tablica wartości górnej przekątnej macierzy 
    //  b[] - tablica wyrazow wolnych b
    //  x[] - tablica rozwiazan

    //  Zwraca: Nic -> operuje na wskaźnikach
    //-------------------------------------------------------------------
    
    
    // Eliminacja w przód dla wektora b; 
    // modyfikujemy tablicę b "w miejscu" - bez alokowania nowej tablicy
    for (int i = 1; i < N; i++) {
        long double m = l[i] / d[i - 1];
        b[i] = b[i] - m * b[i - 1];
    }

    // Podstawianie wsteczne
    x[N - 1] = b[N - 1] / d[N - 1];
    for (int i = N - 2; i >= 0; i--) {
        x[i] = (b[i] - u[i] * x[i + 1]) / d[i];
    }
}



void thomasldpack::Thomas(int N, const long double l[], long double d[], const long double u[],
        long double b[], long double x[]) {
    //-------------------------------------------------------------------
    // Główna funkcja rozwiązująca układ Ax = b dla macierzy trójdiagonalnej
    // l[1..N-1], d[0..N-1], u[0..N-2], b[0..N-1]; wynik w x[0..N-1]
    // UWAGA: d[] i b[] są modyfikowane w miejscu!

    //  Argumenty:
    //  N - rozmiar macierzy A
    //  l[] - tablica wartości dolnej przekątnej macierzy 
    //  d[] - tablica wartości głównej przekątnej macierzy 
    //  u[] - tablica wartości górnej przekątnej macierzy 
    //  b[] - tablica wyrazow wolnych b
    //  x[] - tablica rozwiazan

    //  Zwraca: Nic -> operuje na wskaźnikach
    //-------------------------------------------------------------------

    thomas_procedure_1(N, l, d, u);
    thomas_procedure_2(N, l, u, d, b, x);
}