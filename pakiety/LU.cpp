#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "LU.h"

using namespace std;


void lupack::swap(int* a, int* b) {
//---------------------------------------------------------------------
//  funkcja zamienia ze sobą dwie przekazane wartości
//
//  Argumenty:
//      a - wskaźnik do 1. elementu
//      b - wskaźnik do 2. elementu
//
//  Zwraca:
//      Nic -> funkcja zamienia wartości bezpośrednio w przekazanych elem.
//---------------------------------------------------------------------

    int temp = *a; // Zapisz wartość a w temp
    *a = *b;       // Przypisz wartość b do a
    *b = temp;     // Przypisz temp do b
}


void lupack::LU_decompose(long double A[], int index[], int N){
    //---------------------------------------------------------------------
//  Funkcja dokonująca dekompozycji LU macierzy A z częściowym wyborem
//  elementu podstawowego, ale bez fizycznej zamiany wierszy.
//  
//  Argumenty:
//      A[]             - JEDNOWYMIAROWA tablica zawierająca macierz (porządek wierszowy)
//                        na której wykonamy modyfikacje dekompozycji LU.
//
//      N               - Rozmiar macierzy (tutaj n = N)
//
//      index[]         - Tablica indeksów; początkowo index[i]=i, a później wskaże kolejność
//                        wirtualnych wierszy po zamianach.
//
//  Zwraca: 
//      -Nic -> funkcja zamienia wartości bezpośrednio w przekazanych elem.
//
//
//  ___________________________________________________________________________
//  UWAGA, WAŻNE! DLA OPTYMALIZACJI OBLICZEŃ WSZYSTKIE PRZELICZONE WARTOŚCI SĄ
//  PRZECHOWYWANE W JEDNEJ TABLICY (W PRZYPADKU MACIERZY L JEDYNKI NIEJAWNIE):
//
//  Po wykonaniu funkcji, macierz A zawiera macierze: górnotrójkątną (U) 
//  oraz dolnotrójkątną (mnożniki L, przyjmujemy założenie L[i,i]=1) 
//  ___________________________________________________________________________
//
//  elementy macierzy – pamiętamy, że odczytujemy
//  je zawsze poprzez indeksację: A[index[r]*n+c].
//---------------------------------------------------------------------

    
    
    //  Inicjalizacja tablicy index – początkowy porządek naturalny: 0, 1, 2, ..., n-1
    for (int i = 0; i < N; i++) {
        index[i] = i;
    }
    
    //  Przechodzimy kolumna po kolumnie
    //  k - indeks kolumny na jakiej wykonujemy operacje
    for (int k = 0; k < N; k++) {
        
        //  zamiana wiersza (wybór elementu podstawowego) wykonujemy TYLKO wtedy,
        //  gdy bieżący element podstawowy (A[index[k]*n+k]) jest równy 0.
        if (fabsl(A[index[k] * N + k]) == 0) {
            
            //  Gdy obecny el. podstawowy = 0, to szukany jest inny, największy 
            //  spośród pozostałych elementów w kolumnie k:
            
            int swapIndex = k;
            long double maxVal = fabsl(A[index[k] * N + k]);
            
            for (int i = k + 1; i < N; i++) {
                long double val = fabsl(A[index[i] * N + k]);
                if (val > maxVal) {
                    maxVal = val;
                    swapIndex = i;
                }
            }

            if (maxVal == 0.0L) {printf("\nLU-macierz jest osobilwa/bliska osobilwosci.\n"); exit(1);}
            // Obsłużenie sytuacji, gdy macierz jest osobliwa lub bardzo bliska osobliwości
        

            if (swapIndex != k) {
                lupack::swap(&index[k], &index[swapIndex]);
            }
            //  zamiana kolejności wierszy przy pomocy wektora indeksów
        }
        
        // Eliminacja Gaussa: dla każdego wirtualnego wiersza poniżej (i = k+1,..., n-1)
        for (int i = k + 1; i < N; i++) {
            long double multiplier = A[index[i] * N + k] / A[index[k] * N + k];  
            
            // Zapisujemy mnożnik w miejscu elementu A (on stanowi element L, przyjmujemy L[i,k]=multiplier)
            A[index[i] * N + k] = multiplier;
            // WAŻNE!!!: tam gdzie w eliminacji Gaussa powstałoby 0, zapisujemy mnożnik by nie marnować miejsca
            
            // Aktualizacja pozostałych elementów wirtualnego wiersza
            for (int j = k + 1; j < N; j++) {
                A[index[i] * N + j] -= multiplier * A[index[k] * N + j];
            }
        }
    }
}



void lupack::LU_solve(long double A[], int index[], long double b[], int N) {
//---------------------------------------------------------------------
//  Funkcja rozwiązująca układ równań A*x = b, wykorzystując wcześniej
//  wykonaną dekompozycję LU, przy czym wszystkie operacje odbywają się
//  na "wirtualnych" wierszach określonych przez tablicę index.
//  Proces dzieli się na dwa etapy:
//      1) Podstawienie w przód: rozwiązuje L*y = P*b,
//      2) Podstawienie wstecz: rozwiązuje U*x = y.
//  Wynik (rozwiązanie x) jest zapisywany w tablicy b.
//
//  Argumenty:
//      A[]             - JEDNOWYMIAROWA tablica zawierająca macierz (porządek wierszowy)
//                        na której wykonamy modyfikacje dekompozycji LU.
//      
//      b[]             - tablica wyrazów wolnych (do niej także zapisywane rozwiązania)
//
//      index[]         - Tablica indeksów; początkowo index[i]=i, a później wskaże kolejność
//                        wirtualnych wierszy po zamianach.
//
//      n               - Rozmiar macierzy (tutaj n = N)
//
//  Zwraca: 
//      Nic -> funkcja zamienia wartości bezpośrednio w przekazanych elem.
//---------------------------------------------------------------------
    
    
    // 1. Forward substitution: L * y = P*b
    for (int i = 0; i < N; i++) {
        long double sum = b[index[i]];

        // Odejmujemy wpływ poprzednich elementów (L jest jednostkowa na przekątnej, 
        // a mnożniki są zapisane w A w pozycji: A[index[i] * N + j])
        for (int j = 0; j < i; j++) {
            sum -= A[index[i] * N + j] * b[index[j]];
        }

        // Wynik forward eliminacji zapisujemy z powrotem w wektorze b
        b[index[i]] = sum;
    }


    // 2. Backward substitution: U * x = y
    // Przechodzimy od ostatniego wirtualnego wiersza do pierwszego
    for (int i = N - 1; i >= 0; i--) {
        long double sum = b[index[i]];

        for (int j = i + 1; j < N; j++) {
            sum -= A[index[i] * N + j] * b[index[j]];
        }

        // Uzyskujemy dzieląc przez przekątny element macierzy U (przechowywany w A)
        sum /= A[index[i] * N + i];
        b[index[i]] = sum;
    }

}



void lupack::LU_decompose_and_solve(long double A[], long double b[], int N){
//---------------------------------------------------------------------
//  Funkcja dekomponuje przekazaną macierz na macierze L oraz U,
//  a następnie rozwiązuje układ równań
//
//  Argumenty:
//      A[]             - JEDNOWYMIAROWA tablica zawierająca macierz (porządek wierszowy)
//      
//      b[]             - tablica wyrazów wolnych (do niej także zapisywane rozwiązania)
//
//      N              - Rozmiar macierzy (tutaj n = N)
//
//  Zwraca: 
//      Nic -> funkcja zamienia wartości bezpośrednio w przekazanych elem.
//---------------------------------------------------------------------
    
    int* index = new int[N];

    lupack::LU_decompose(A,index, N);
    lupack::LU_solve(A, index, b, N);

    delete[] index;
}