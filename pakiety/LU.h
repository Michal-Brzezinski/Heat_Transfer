#ifndef __lu_h
#define __lu_h

namespace lupack{

    void swap(int* a, int* b);
    void LU_decompose(long double A[], int index[], int n);
    void LU_solve(long double A[], int index[], long double b[], int n);
    void LU_decompose_and_solve(long double A[], long double b[], int n);
}

#endif