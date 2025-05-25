#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "LU.h"

using namespace std;


void ludecomposepack::swap(int* a, int* b) {
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


void ludecomposepack::LU_decompose(long double A[], int index[], int n){
    //---------------------------------------------------------------------
//  Funkcja dokonująca dekompozycji LU macierzy A z częściowym wyborem
//  elementu podstawowego, ale bez fizycznej zamiany wierszy.
//  
//  Argumenty:
//      A[]             - JEDNOWYMIAROWA tablica zawierająca macierz (porządek wierszowy)
//                        na której wykonamy modyfikacje dekompozycji LU.
//
//      n               - Rozmiar macierzy (tutaj n = N)
//
//      index[]         - Tablica indeksów; początkowo index[i]=i, a później wskaże kolejność
//                        wirtualnych wierszy po zamianach.
//
//  Zwraca: 
//      -1 lub 0 w zależności od powodzenia operacji
//
//
//  UWAGA WAŻNE!!!!! DLA OPTYMALIZACJI OBLICZEŃ WSZYSTKIE PRZELICZONE WARTOŚCI SĄ
//  PRZECHOWYWANE W JEDNEJ TABLICY (W PRZYPADKU MACIERZY L JEDYNKI NIEJAWNIE):
//
//  Po wykonaniu funkcji, macierz A zawiera górną (U) oraz dolną (mnożniki L, 
//  przyjmujemy założenie L[i,i]=1) 
//
//  elementy macierzy – pamiętamy, że odczytujemy
//  je zawsze poprzez indeksację: A[index[r]*n+c].
//---------------------------------------------------------------------

    
    
    // Inicjalizacja tablicy index – początkowy porządek naturalny: 0, 1, 2, ..., n-1
    for (int i = 0; i < n; i++) {
        index[i] = i;
    }
    
    // Przechodzimy kolumna po kolumnie
    // k - indeks kolumny na jakiej wykonujemy operacje
    for (int k = 0; k < n; k++) {
        
        // zamiana wiersza (wybór elementu podstawowego) wykonujemy TYLKO wtedy,
        // gdy bieżący element podstawowy (A[index[k]*n+k]) jest równy 0.
        if (fabsl(A[index[k] * n + k]) == 0) {
            int swapIndex = k;
            long double maxVal = fabsl(A[index[k] * n + k]);
            
            for (int i = k + 1; i < n; i++) {
                long double val = fabsl(A[index[i] * n + k]);
                if (val > maxVal) {
                    maxVal = val;
                    swapIndex = i;
                }
            }

            if (maxVal == 0.0L) exit(1);
            // Macierz jest osobliwa lub bardzo bliska osobliwości,
            // obsłuż błąd, np. zwracając kod błędu lub kończąc działanie.
        

            if (swapIndex != k) {
                ludecomposepack::swap(&index[k], &index[swapIndex]);
            }
            //  zamiana kolejności wierszy przy pomocy wektora indeksów
        }
        
        // Eliminacja Gaussa: dla każdego wirtualnego wiersza poniżej (i = k+1,..., n-1)
        for (int i = k + 1; i < n; i++) {
            long double multiplier = A[index[i] * n + k] / A[index[k] * n + k];  
            // Zapisujemy mnożnik w miejscu elementu A (on stanowi element L, przyjmujemy L[i,k]=multiplier)
            A[index[i] * n + k] = multiplier;
            // WAŻNE!!!: tam gdzie w eliminacji Gaussa powstałoby 0, zapisujemy mnożnik by nie marnować miejsca
            
            // Aktualizacja pozostałych elementów wirtualnego wiersza
            for (int j = k + 1; j < n; j++) {
                A[index[i] * n + j] -= multiplier * A[index[k] * n + j];
            }
        }
    }
}



void ludecomposepack::LU_solve(long double A[], int index[], long double b[], int n) {
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
    for (int i = 0; i < n; i++) {
        long double sum = b[index[i]];
        
        // Odejmujemy wpływ poprzednich elementów (L jest jednostkowa na przekątnej, 
        // a mnożniki są zapisane w A w pozycji: A[index[i]*n+j])
        for (int j = 0; j < i; j++) {
            sum -= A[index[i] * n + j] * b[index[j]];
        }
        b[index[i]] = sum;  // wynik forward eliminacji zapisujemy z powrotem w wektorze b
    }

    // 2. Backward substitution: U * x = y
    // Przechodzimy od ostatniego wirtualnego wiersza do pierwszego
    for (int i = n - 1; i >= 0; i--) {
        long double sum = b[index[i]];
        for (int j = i + 1; j < n; j++) {
            sum -= A[index[i] * n + j] * b[index[j]];
        }
        // uzyskujemy dzieląc przez przekątny element macierzy U (przechowywany w A)
        sum /= A[index[i] * n + i];
        b[index[i]] = sum;
    }
}



void ludecomposepack::LU_decompose_and_solve(long double A[], long double b[], int n){
//---------------------------------------------------------------------
//  Funkcja dekomponuje przekazaną macierz na macierze L oraz U,
//  a następnie rozwiązuje układ równań
//
//  Argumenty:
//      A[]             - JEDNOWYMIAROWA tablica zawierająca macierz (porządek wierszowy)
//      
//      b[]             - tablica wyrazów wolnych (do niej także zapisywane rozwiązania)
//
//      n               - Rozmiar macierzy (tutaj n = N)
//
//  Zwraca: 
//      Nic -> funkcja zamienia wartości bezpośrednio w przekazanych elem.
//---------------------------------------------------------------------
    
    int* index = new int[n];

    ludecomposepack::LU_decompose(A,index, n);
    ludecomposepack::LU_solve(A, index, b, n);

    delete[] index;
}