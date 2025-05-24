#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <utility> // dla std::swap

//----------------------------------------------------------------------
// Stałe fizyczne i przestrzenne
//----------------------------------------------------------------------
const long double D     = 1.0e0L;    // współczynnik dyfuzji
const long double b     = 1.0e-1L;   // parametr początkowy
const long double t_max = 1.0e0L;    // czas końcowy
// a = 6 * sqrt(D*t_max)
const long double a     = 6.0e0L;

//----------------------------------------------------------------------
// Rozmiary siatki
//----------------------------------------------------------------------
//  liczba węzłów siatki przestrzennej
const unsigned long long N = 12000000001ULL;

//  liczba węzłów siatki czasowej 
const unsigned long long M = 10000000001ULL;

//----------------------------------------------------------------------
// Wartości kroków na siatce czasowo-przestrzennej
//----------------------------------------------------------------------
const long double h  = (2*a)/(N-1ULL);     // krok przestrzenny
const long double dt = t_max/(M-1ULL);     // krok czasowy (krok całkowania)



//----------------------------------------------------------------------
// Funkcje pomocnicze
//----------------------------------------------------------------------


void warunek_poczatkowy(long double* U, const long double* X) {
    //-------------------------------------------------------------------
    // Warunek początkowy U(x,0):
    // Funkcja inicjalizuje wartości dla tablicy U na odpowiednie 
    // biorąc pod uwagę podany warunek początkowy
    //-------------------------------------------------------------------

    for (unsigned long long i = 0ULL; i < N; ++i) {
        U[i] = (X[i] < 0.0L) ? 0.0L : std::expl(-X[i] / b);
    }
}


long double analytic_solution(long double x, long double t) {
    //-------------------------------------------------------------------
    // Rozwiązanie analityczne U(x,t)
    //-------------------------------------------------------------------

    long double z = (2.0L * D * t / b - x) / (2.0L * std::sqrtl(D * t));
    long double pref = 0.5L * std::expl(D * t / (b * b) - x / b);
    return pref * std::erfcl(z);
}

// Jawna metoda KMB
void explicit_step(const long double* U_old, long double* U_new, long double lambda) {
    // Warunki brzegowe: U_new[0]=U_new[N-1]=0 (przyjmujemy, że już są ustawione)
    for (unsigned long long i = 1; i + 1 < N; ++i) {
        U_new[i] = U_old[i] + lambda * (U_old[i + 1] - 2.0L * U_old[i] + U_old[i - 1]);
    }
}

// Algorytm Thomasa dla macierzy trójdiagonalnej
void thomas(const long double* aa, const long double* bb, const long double* cc,
            const long double* dd, long double* U) {
    // Alokujemy tymczasowe tablice dla cp i dp
    long double* cp = new long double[N];
    long double* dp = new long double[N];
    
    cp[0] = cc[0] / bb[0];
    dp[0] = dd[0] / bb[0];
    for (unsigned long long i = 1; i < N; ++i) {
        long double m = bb[i] - aa[i] * cp[i - 1];
        cp[i] = ((i + 1 < N) ? cc[i] : 0.0L) / m;
        dp[i] = (dd[i] - aa[i] * dp[i - 1]) / m;
    }
    U[N - 1] = dp[N - 1];
    for (unsigned long long i = N - 1; i-- > 0;) {
        U[i] = dp[i] - cp[i] * U[i + 1];
    }
    
    delete[] cp;
    delete[] dp;
}

// Schemat Laasonena (BTCS)
void laasonen_step(const long double* U_old, long double* U_new, long double lambda) {
    // Alokujemy tablice na współczynniki układu trójdiagonalnego
    long double* aa = new long double[N];
    long double* bb = new long double[N];
    long double* cc = new long double[N];
    long double* dd = new long double[N];
    
    for (unsigned long long i = 0; i < N; ++i) {
        if (i == 0 || i + 1 == N) {
            // Dla warunków Dirichleta U=0
            aa[i] = 0.0L;
            bb[i] = 1.0L;
            cc[i] = 0.0L;
            dd[i] = 0.0L;
        } else {
            aa[i] = -lambda;
            bb[i] = 1.0L + 2.0L * lambda;
            cc[i] = -lambda;
            dd[i] = U_old[i];
        }
    }
    thomas(aa, bb, cc, dd, U_new);
    
    delete[] aa;
    delete[] bb;
    delete[] cc;
    delete[] dd;
}

// Maksymalny błąd między rozwiązaniem numerycznym a analitycznym w danym czasie t
long double compute_max_error(const long double* U_num, const long double* X, long double t) {
    long double max_err = 0.0L;
    for (unsigned long long i = 0; i < N; ++i) {
        long double ue = analytic_solution(X[i], t);
        long double e = std::fabsl(U_num[i] - ue);
        if (e > max_err) {
            max_err = e;
        }
    }
    return max_err;
}

int main() {
    // Wypisanie wymiarów siatki
    std::cout << "N = " << N << ", M = " << M << std::endl;

    // Alokacja tablic dynamicznych
    long double* X   = new long double[N];
    long double* U0  = new long double[N];
    long double* Ue  = new long double[N];
    long double* Ul  = new long double[N];
    long double* Tmp = new long double[N];

    // Utworzenie siatki przestrzennej jako: X[i] = -a + i*h
    for (unsigned long long i = 0; i < N; ++i) {
        X[i] = -a + static_cast<long double>(i) * h;
    }
    
    // Inicjalizacja warunku początkowego U(x,0)
    warunek_poczatkowy(U0, X);
    // Kopiujemy U0 do Ue oraz Ul
    for (unsigned long long i = 0ULL; i < N; ++i) {
        Ue[i] = U0[i];
        Ul[i] = U0[i];
    }


    //  definiujemy parametr lambda dla obu metod
    long double lambda = D * dt / (h * h);

    // Pętla czasowa
    for (unsigned long long n = 0ULL; n < M; ++n) {
        // Metoda KMB
        explicit_step(Ue, Tmp, lambda);
        // Zamiana wskaźników, aby uniknąć kopiowania tablic – teraz Ue wskazuje na wynik nowej iteracji
        std::swap(Ue, Tmp);
        
        // Schemat Laasonena (implicit)
        laasonen_step(Ul, Tmp, lambda);
        std::swap(Ul, Tmp);
    }

    // Obliczenie błędów przy czasie t_max
    long double err_e = compute_max_error(Ue, X, t_max);
    long double err_l = compute_max_error(Ul, X, t_max);
    std::cout << "Max error explicit = " << err_e << std::endl;
    std::cout << "Max error laasonen = " << err_l << std::endl;

    // Zapis wyników do pliku CSV
    std::ofstream fout("results.csv");
    fout << "x,U_explicit,U_laasonen,U_exact\n";
    for (unsigned long long i = 0; i < N; ++i) {
        long double ue_x = analytic_solution(X[i], t_max);
        fout << X[i] << "," << Ue[i] << "," << Ul[i] << "," << ue_x << "\n";
    }
    fout.close();

    // Dealokacja pamięci
    delete[] X;
    delete[] U0;
    delete[] Ue;
    delete[] Ul;
    delete[] Tmp;

    return 0;
}
