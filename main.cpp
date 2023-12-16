#include <iostream>
#include <cmath>
using namespace std;

#include "struktury.h"
#include "ElementUniversal.h"

void wczytywanie(Grid& siatka);

double fun1(double x){  double w =5*pow(x,2)+3*x+6;    return w;  }
double fun2(double x, double y){  double w =5*pow(x,2)*pow(y,2)+3*x*y+6;    return w;  }


int main() {
   Grid siatka;

   wczytywanie(siatka);
    cout<<endl;
   // wyswietl_nodes(siatka);
  // Gauss(4).wym1D(fun1);  // Gauss(4).wym2D(fun2);

   int ile_pkt=2;
   ElementUniversal elementUniversal(ile_pkt);
   elementUniversal.obliczanie_ksi_etta();
   elementUniversal.wyswietl();


    cout<<endl;
    elementUniversal.oblicanie_tabN();
    for (int i = 0; i < 4; ++i) {
        cout<<"surface "<<i+1<<endl;
        for (int j = 0; j < ile_pkt; ++j) {
            for (int k = 0; k < 4; ++k) {
                cout<<setw(15)<<elementUniversal.surface[i].N[j][k]<<"  ";
            }cout<<endl;
        }cout<<endl;
    }

    elementUniversal.obliczanie_tab_FN();
    elementUniversal.wyswietl_tab_FN();


    for (int i = 0; i < dane.Elements_number; ++i) {
       cout<<"****************************************************************************"<<endl;
       cout<<"Element "<<i+1<<endl;
       cout<<"****************************************************************************"<<endl;
       elementUniversal.jakobian(siatka.Elements[i], dane.Conductivity,i);cout<<endl;
       elementUniversal.macierz_Hbc(ile_pkt,siatka.Elements[i],i);//i wektor P
        siatka.Elements[i].obliczanie_h_calk();
        cout<<endl;
        cout<<"-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --"<<endl;
        cout<<"H calkowite elmentu "<<i+1<<endl;
        siatka.Elements[i].display_H_calk();  cout << endl;
        cout<<"P calkowite elmentu "<<i+1<<endl;
        siatka.Elements[i].display_P_calk();  cout << endl;
        cout<<"C calkowite elmentu "<<i+1<<endl;
        siatka.Elements[i].display_C_calk();
        cout<<"-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --"<<endl;
   }

    SOE soe;
    soe.agregacja(siatka);
    cout<<endl<<"Globalne H"<<endl;soe.display_G_H();
    cout<<endl<<"Globalne P"<<endl;soe.display_G_P();
    cout<<endl<<"Globalne C"<<endl;soe.display_G_C();

    soe.rozwiazywanie_temp();


/*  //wszytskie naraz h calk
    for (int i = 0; i < dane.Elements_number; ++i) {
        cout<<"H calkowite elmentu "<<i<<endl;
        siatka.Elements[i].display_H_calk();
    }*/


/* cout<<"\n\n****************************************************************************\n"

       "Sprawdzenie czy h i hbc sa w struct element"<<endl;
        for (int q = 0; q < dane.Elements_number; ++q) {
            cout << "element " << q + 1 << endl<<endl;
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    cout << siatka.Elements[q].tab_H[i][j] << "   ";
                }
                cout << endl;
            }cout << endl;
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    cout << siatka.Elements[q].Hbc[i][j] << "   ";
                }
                cout << endl;
            } cout << endl;

        }*/



  return 0;
}


