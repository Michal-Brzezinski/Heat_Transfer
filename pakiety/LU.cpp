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
//  Funkcja rozwiązująca układ równań Ax = b, wykorzystując wcześniej
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
// Ten fragment kodu rozwiązuje układ równań L*y = P*b. 
// Tutaj L jest macierzą dolnotrójkątną z jedynkami na przekątnej (nie zapisujemy jej elementów, 
// gdyż przekątna jest domyślnie równa 1), a P*b oznacza wektor b przestawiony zgodnie z
// permutacjami ustalonymi podczas dekompozycji LU. Wektor 'index' przechowuje numery 
// wierszy odpowiadające nowej kolejności, dzięki czemu mamy dostęp do prawidłowych wierszy macierzy i wektora.
for (int i = 0; i < N; i++) {
    // Inicjujemy zmienną 'sum' wartością elementu b odpowiadającą bieżącemu równaniu, 
    // już przestawionego przez pivoting (indeksowanym przez index[i]).
    long double sum = b[index[i]];

    // Odejmujemy wpływ poprzednich elementów (L jest jednostkowa na przekątnej, 
    // a mnożniki są zapisane w A w pozycji: A[index[i] * N + j])
    // Następnie dla każdego wcześniejszego równania (dla j od 0 do i-1)
    // odejmujemy iloczyn odpowiadającego współczynnika z macierzy L oraz już wyznaczonej 
    // niewiadomej y (przechowywanej w b[index[j]]).
    // Dzięki temu eliminujemy wpływ wcześniej obliczonych wartości na bieżące równanie, 
    // uzyskując w efekcie wartość y dla wiersza index[i].
    for (int j = 0; j < i; j++) {
        sum -= A[index[i] * N + j] * b[index[j]];
    }

    // Wynik forward eliminacji zapisujemy z powrotem w wektorze b
    // Przypisanie sum do b[index[i]] oznacza, że obliczone rozwiązanie 
    // dla y danego równania zastępuje oryginalną wartość z wektora b.
    // W kolejnych operacjach b będzie już zawierało elementy y.
    b[index[i]] = sum;
}


// 2. Backward substitution: U * x = y
// Przechodzimy od ostatniego wirtualnego wiersza do pierwszego
// Ten fragment kodu rozwiązuje układ równań U*x = y, gdzie U jest macierzą górnotrójkątną 
// uzyskaną z dekompozycji LU. Algorytm zaczyna od ostatniego równania (które zawiera tylko jedno niewiadome)
// i idzie "w górę", eliminując wpływ już wyznaczonych zmiennych, aby uzyskać ostateczne wartości x.
// Wektor 'index' jest ponownie używany, by odzwierciedlić kolejność wierszy.
for (int i = N - 1; i >= 0; i--) {
    // Inicjujemy zmienną 'sum' wartością odpowiadającą elementowi y, który wcześniej zapisaliśmy 
    // w wektorze b w wyniku forward substitution.
    long double sum = b[index[i]];

    // Iterujemy po kolumnach dla każdego równania (od i+1 do N-1), 
    // czyli po elementach macierzy U należących do prawej strony przekątnej.
    // Odejmujemy iloczyn współczynnika U z wcześniej obliczoną 
    // niewiadomą x (przechowywaną w b[index[j]]) od skumulowanej sumy.
    // W ten sposób eliminujemy znane już składniki, pozostawiając tylko wyraz związany z x[i].
    for (int j = i + 1; j < N; j++) {
        sum -= A[index[i] * N + j] * b[index[j]];
    }

    // Uzyskujemy dzieląc przez przekątny element macierzy U (przechowywany w A)
    // Każde równanie układu U*x = y ma postać:
    // A[index[i]*N + i] * x[index[i]] + (suma iloczynów z już znanymi x) = y[index[i]]
    // Dlatego po wyeliminowaniu pozostałych składników dzielimy skorygowaną 
    // sum przez element diagonalny, aby uzyskać wartość x dla bieżącego równania.
    sum /= A[index[i] * N + i];
    
    // Wynik rozwiązania dla x zapisujemy z powrotem w wektorze b, 
    // nadpisując tym samym wcześniej obliczone y. W efekcie, po zakończeniu pętli,
    // wektor b zawiera ostateczne rozwiązanie układu równań A*x = b.
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