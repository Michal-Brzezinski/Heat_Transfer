#ifndef __thomas_h
#define __thomas_h

namespace thomaspack{

    void thomas_procedure_1(int N, const long double l[], long double d[], const long double u[]);
    
    void thomas_procedure_2(int N, const long double l[], const long double u[], 
        const long double d[], long double b[], long double x[]);
    
    void Thomas(int N, const long double l[], long double d[], 
        const long double u[], long double b[], long double x[]);

}

#endif