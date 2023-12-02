#ifndef MES_STRUKTURY_H
#define MES_STRUKTURY_H
#include <iostream>
#include <iomanip>
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

    Element(){ ide=0; for (int i = 0; i < 4; i++) { ID[i]=Node();}   }

    void set_tabH(double tab[4][4]);
    void set_Hbc(double tab[4][4]);
    void set_P(double tab[4]);
    void display_H_calk();

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
    SOE();
    void agregacja(Grid &siatka);
    void display_G_H();
    void display_G_P();

  //  ~SOE();
};

void wyswietl_dane();
void wyswietl_nodes(Grid siatka);
void wyswietl_elements(Grid siatka);
void wyswietl_BC(Grid siatka);

#endif //MES_STRUKTURY_H
