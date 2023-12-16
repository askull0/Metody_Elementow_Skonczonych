#ifndef MES_ELEMENTUNIVERSAL_H
#define MES_ELEMENTUNIVERSAL_H
#include "struktury.h"
#include "Gauss.h"

class Surface {
public:
    int e;
    double **N;
   Surface(){};//Surface():e(0), N(nullptr){ };
   Surface(int e);

    //~Surface();
};


class ElementUniversal{
public:
    int E;
    Gauss obj;
    double **dtksi;//macierz przechwujaca wartossci po etta
    double **dtetta;//macierz przechwujaca wartossci po ksi

    double **tab_FN;//wartosci f. ksztaltu w pkt calk

    double H_kon[4][4];
    double C_kon[4][4];
    double Hbc_kon[4][4];
    Surface surface[4];

    ElementUniversal(int e);
//dNksi - e , przesylane etta- n
    double dNksi1(double n);
    double dNksi2(double n);
    double dNksi3(double n);
    double dNksi4(double n);

    double dNetta1(double e);
    double dNetta2(double e);
    double dNetta3(double e);
    double dNetta4(double e);

    void obliczanie_ksi_etta();
    void wyswietl();

    void obliczanie_tab_FN();
    void wyswietl_tab_FN();

    void jakobian(Element& elem, int& GD_k,int nr_e);  //liczy H I C

    void display_tab_22(double tab[2][2]);
    void display_tab_4(double tab[4]);
    void display_tab_44(double tab[4][4]);
    void display_tab_E4(double **tab);
    void display_H(double tab[4][4]);

    //surface
    void oblicanie_tabN();
    void funkcja_do_tabN(int i , double a, double b, int q);

    void macierz_Hbc(int x,Element& elem,int nr_e);   //liczy tez wektor P

    double F_N1(double ksi, double etta);
    double F_N2(double ksi, double etta);
    double F_N3(double ksi, double etta);
    double F_N4(double ksi, double etta);

    ~ElementUniversal();

};



#endif //MES_ELEMENTUNIVERSAL_H
