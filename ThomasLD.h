#ifndef __thomasld_h
#define __thomasld_h

namespace thomasldpack{

    void thomas_procedure_1(unsigned long long N, const long double l[], long double d[], const long double u[]);
    
    void thomas_procedure_2(unsigned long long N, const long double l[], const long double u[], 
        const long double d[], long double b[], long double x[]);
    
    void Thomas(unsigned long long N, const long double l[], long double d[], 
        const long double u[], long double b[], long double x[]);

}

#endif