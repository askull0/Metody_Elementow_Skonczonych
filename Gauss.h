#ifndef MES_GAUSS_H
#define MES_GAUSS_H
#include <cmath>
#include <iostream>
using namespace std;

class Gauss{
public:
    int n; // wybor pkt
    double *w, *x;

    Gauss (int N);
    void wym1D(double(*f)(double));
    void wym2D(double(*f)(double,double));
    ~Gauss();
};


#endif //MES_GAUSS_H
