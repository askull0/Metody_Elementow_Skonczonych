#ifndef MES_STRUKTURY_H
#define MES_STRUKTURY_H
#include <iostream>
#include <iomanip>
#include <cmath>
const double eps = 1e-12;
using namespace std;

struct GlobalData{
    int SimulationTime, SimulationStepTime, Conductivity, Alfa, Tot,
            InitialTemp, Density, SpecificHeat, Nodes_number, Elements_number;
};
extern GlobalData dane;

struct Node{
    int id;
    double X,Y;
    int BC;
    Node(){id=0.;X=0.;Y=0.;BC=0;}
};

class Element {
public:
    int ide;
    Node ID[4];
    double tab_H[4][4];
    double Hbc[4][4];
    double P[4];
    double H_CALK[4][4];  // suma H i Hbc
    double C[4][4];
    double t0[4];  //wektor t0 z temp poczatkową
    double C_tau[4][4];  //macierz C podzielona przez tau*/
    double C_t0[4];   //C_tau razy t0


    Element(){ ide=0; for (int i = 0; i < 4; i++) { ID[i]=Node();}  // fill_n(&C_tau[0][0], 4 * 4, 0);
          for (int i = 0; i < 4; ++i) {  t0[i]=dane.InitialTemp;  }
    }

    void set_tabH(double tab[4][4]);
    void set_Hbc(double tab[4][4]);
    void set_P(double tab[4]);
    void set_C(double tab[4][4]);
    void display_H_calk();
    void display_P_calk();

    void obliczanie_h_calk();

};

class Grid {
public:
    Node* Nodes;
    Element* Elements;

    Grid();
   // ~Grid();
};

class SOE{
public:
    double **G_H;
    double *G_P;
    double G_HP;  //chyba test do innego rozw ukladu rownan

    SOE();
    void agregacja(Grid &siatka);
    void display_G_H();
    void display_G_P();

    //uklad rownan

    static void print(double ** A);
    static void swap_row(double **A, double * B, int i, int j);
    static int eliminacja(double **A,  double *B);
    static double* result(double **A,  double *B, double *x);
    static void gauss_elim(double **A,  double *B, double *x);
    void gauss_crout(double *x);

    bool gauss_crout2(double *X);

  //  ~SOE();
};

void wyswietl_dane();
void wyswietl_nodes(Grid siatka);
void wyswietl_elements(Grid siatka);
void wyswietl_BC(Grid siatka);

#endif //MES_STRUKTURY_H
