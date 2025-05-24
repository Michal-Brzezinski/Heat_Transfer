#ifndef __thomasld_h
#define __thomasld_h

namespace thomasldpack{
    
    void thomas_procedure_1(int N, const double l[], double d[], const double u[]);
    
    void thomas_procedure_2(int N, const double l[], const double u[], 
        const double d[], double b[], double x[]);
    
        void Thomas(int N, const double l[], double d[], const double u[],
        double b[], double x[]);

}

#endif