#include "ElementUniversal.h"

Surface::Surface(int e) :e(e){
    N = new double*[e];
    for (int i = 0; i < e; ++i) {
        N[i] = new double[4];
        for (int j = 0; j < 4; ++j) {
            N[i][j] = 0.0;
        }
    }
}

ElementUniversal::ElementUniversal(int e) : E(pow(e,2)), obj(e) { //bedzie liczyc dla jakiego elemntu calke
    dtetta=new double*[E];
    dtksi=new double*[E];
    tab_FN=new double*[E];

    for (int i = 0; i < 4; ++i) {
        surface[i] = Surface( e );
    }

    for (int i = 0; i < E; i++) {
        dtetta[i] = new double[4];
        dtksi[i] = new double[4];// tablica ksi przeyslamy etta-n, y // dN/dksi
        tab_FN[i] = new double[4];
    }

}//dla jakiego elemntu calke czyli np. 2 pkt lub 3 lub 4

//dNksi - e , przesylane etta- n
double  ElementUniversal::dNksi1(double n){return (-1./4.*(1-n));}
double  ElementUniversal::dNksi2(double n){return (-1./4.*(1+n));}
double  ElementUniversal::dNksi3(double n){return (1./4.*(1+n));}
double  ElementUniversal::dNksi4(double n){return (1./4.*(1-n));}

double  ElementUniversal::dNetta1(double e){return (-1./4.)*(1-e); }
double  ElementUniversal::dNetta2(double e){return (1./4.)*(1-e);}
double  ElementUniversal::dNetta3(double e){return (1./4.)*(1+e);}
double  ElementUniversal::dNetta4(double e){return (-1./4.)*(1+e);}

void  ElementUniversal::obliczanie_ksi_etta(){//oblczianie tabelek z wartosciami po ksi i etta
    for (int i = 0; i < 4; i++) {
        int licznik = 0, licznik2 = 0;
        for (int j = 0; j < E; j++) {
            if (licznik <= obj.n) {
                if (licznik == obj.n) {
                    licznik = 0;
                    licznik2++;
                }
                switch (i) {//blokownie kolumn, liczenie najpierw tylko po funkacj kszt czyli i=0 dN1/dksi
                    case 0:
                        dtksi[j][i] = dNetta1(obj.x[licznik2]);
                        dtetta[j][i] = dNksi1(obj.x[licznik]);
                        break;
                    case 1:
                        dtksi[j][i] = dNetta2(obj.x[licznik2]);
                        dtetta[j][i] = dNksi2(obj.x[licznik]);
                        break;
                    case 2:
                        dtksi[j][i] = dNetta3(obj.x[licznik2]);
                        dtetta[j][i] = dNksi3(obj.x[licznik]);
                        break;
                    case 3:
                        dtksi[j][i] = dNetta4(obj.x[licznik2]);
                        dtetta[j][i] = dNksi4(obj.x[licznik]);
                        break;
                }
                licznik++;
            }
        }
    }
}

void  ElementUniversal::wyswietl(){
    cout<<"****************************************************************************"<<endl;
    cout<<"Tablica po dN/detta"<<endl;  //sprawdz
    display_tab_E4(dtetta);
    cout<<endl<<"Tablica po dN/dksi"<<endl;
    display_tab_E4(dtksi);
}

void ElementUniversal::obliczanie_tab_FN(){
    for (int i = 0; i < 4; ++i) {
        int a=0,b=0; //liczniki
        for (int j = 0; j < E; ++j) {
            if (a <= obj.n) {
                if (a == obj.n) {  a = 0;  b++;   }
                switch(i){
                    case 0:
                        tab_FN[j][i]=F_N1(obj.x[a], obj.x[b]);break;
                    case 1:
                        tab_FN[j][i]=F_N2(obj.x[a], obj.x[b]);break;
                    case 2:
                        tab_FN[j][i]=F_N3(obj.x[a], obj.x[b]);break;
                    case 3:
                        tab_FN[j][i]=F_N4(obj.x[a], obj.x[b]);break;
                }
                a++;
            }
        }
    }

}
void ElementUniversal::wyswietl_tab_FN(){
    cout<<"****************************************************************************"<<endl;
    cout<<"wartosci FN"<<endl;
    display_tab_E4(tab_FN);
}

void ElementUniversal::jakobian(Element& elem, int& GD_k,int nr_e){  //macierz H i macierz C

    double tab_jakobian[2][2], tab_j_wyzn[2][2];//jakobian przemnozony przez odwrotnosc wyzn_J
    double **tab_X, **tab_Y;
    double tmp[4][4], tmpx[4][4], tmpy[4][4];
    double H_suma[4][4], H_pkt[4][4];
    int count_n=0,count_m=0;//liczniki
    double wyz_j=0.;
    double tab_C[4][4];    fill_n(&tab_C[0][0], 4 * 4, 0);
    double C_pkt[4][4];    fill_n(&C_pkt[0][0], 4 * 4, 0);
    double tmp_C[4][4];      fill_n(&tmp[0][0], 4 * 4, 0);

    tab_X = new double*[E];//ile wierszy
    tab_Y = new double*[E];//wierszy

    for (int i = 0; i < E; i++) {
        tab_X[i] = new double[4];
        tab_Y[i] = new double[4];
    }

    for(int q=0; q<E; q++) { // q - odp konrektnemu pkt calkowanie

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                for (int k = 0; k < 4; k++) {
                    if(i==0 & j==0){ tab_jakobian[0][0] += (dtksi[q][k] *  elem.ID[k].X);}  //       punkty[k].X
                    if(i==0 & j==1){ tab_jakobian[0][1] += (dtksi[q][k] *  elem.ID[k].Y);}  //       punkty[k].Y
                    if(i==1 & j==0){ tab_jakobian[1][0] += (dtetta[q][k] * elem.ID[k].X);}  //       punkty[k].X
                    if(i==1 & j==1){ tab_jakobian[1][1] += (dtetta[q][k] * elem.ID[k].Y);}  //       punkty[k].Y
                }
            }
        }
    //    cout<<endl<<"Jakobian dla "<<q+1<<" pkt calk"<<endl;
    //    display_tab_22(tab_jakobian);
        wyz_j = tab_jakobian[0][0] * tab_jakobian[1][1] - tab_jakobian[0][1] * tab_jakobian[1][0];
    //    cout<<"det J = "<<wyz_j<<endl;

        //jakobian przemnozony i odwrocony
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                if(i==0 & j==0){tab_j_wyzn[0][0] = (1./wyz_j)*tab_jakobian[1][1]; }
                if(i==0 & j==1){tab_j_wyzn[0][1] = (1./wyz_j)*(-1.)*tab_jakobian[0][1]; }
                if(i==1 & j==0){tab_j_wyzn[1][0] = (1./wyz_j)*(-1.)*tab_jakobian[1][0];   }
                if(i==1 & j==1){tab_j_wyzn[1][1] = (1./wyz_j)*tab_jakobian[0][0];    }
            }
        }
        // cout<<endl<<"Jakobian przemnozony"<<endl;display_tab_22(tab_j_wyzn);

        //dN/x i dN/y
        for (int j = 0; j < 4; ++j) {
            tab_X[q][j]=tab_j_wyzn[0][0]*dtksi[q][j]+tab_j_wyzn[0][1]*dtetta[q][j];
            tab_Y[q][j]=tab_j_wyzn[1][0]*dtksi[q][j]+tab_j_wyzn[1][1]*dtetta[q][j];
        }

        //macierz razy m.transponowane i wyliczanie H dla pkt  oraz wyliczanie C dla pkt
        if (count_m <= sqrt(E)) {
            if (count_m == sqrt(E)) {
                count_m = 0;
                count_n++;
            }
        }

        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                tmpx[i][j] = tab_X[q][i] * tab_X[q][j];
                tmpy[i][j] = tab_Y[q][i] * tab_Y[q][j];
                H_suma[i][j]=tmpx[i][j]+tmpy[i][j];
                H_pkt[i][j]=GD_k*(H_suma[i][j])*(wyz_j);

                tmp[i][j] = H_pkt[i][j]*obj.w[count_m] * obj.w[count_n];
                H_kon[i][j]+=tmp[i][j];

                //macierz C
                tmp_C[i][j] = tab_FN[q][i] * tab_FN[q][j];
                C_pkt[i][j]=tmp_C[i][j]*dane.Density*dane.SpecificHeat*wyz_j;
                tab_C[i][j] = C_pkt[i][j]*obj.w[count_m] * obj.w[count_n];

                C_kon[i][j]+=tab_C[i][j];

            }
        } count_m++;

        //wyczyszczenie tablicy
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                tab_jakobian[i][j]=0;
            }
        }

    }//koniec petli q dla ilosci pkt c

    cout<<endl << "tab_X - po dN/x" << endl;
    display_tab_E4(tab_X);
    cout << "tab_Y - po dN/y" << endl;
    display_tab_E4(tab_Y);

    cout<<endl<<"H dla elementu "<<nr_e+1<<endl;
    display_H(H_kon);
    elem.set_tabH(H_kon);

    cout<<endl<<"C dla elementu "<<nr_e+1<<endl;
    display_tab_44(C_kon);
    elem.set_C(C_kon);

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            C_kon[i][j]=0.;
            H_kon[i][j]=0.;
        }
    }
}


double  ElementUniversal::F_N1(double ksi, double etta){ return(1./4.*(1.-ksi)*(1.-etta)); }
double  ElementUniversal::F_N2(double ksi, double etta){ return(1./4.*(1.+ksi)*(1.-etta)); }
double  ElementUniversal::F_N3(double ksi, double etta){ return(1./4.*(1.+ksi)*(1.+etta)); }
double  ElementUniversal::F_N4(double ksi, double etta){ return(1./4.*(1.-ksi)*(1.+etta)); }

void ElementUniversal::funkcja_do_tabN(int i ,double a, double b, int q) {
        for (int j = 0; j < 4; ++j) {
            switch (j) {
                case 0:
                    surface[q].N[i][j]=F_N1(a, b);break;
                case 1:
                    surface[q].N[i][j]=F_N2(a,b);break;
                case 2:
                    surface[q].N[i][j]=F_N3(a,b);break;
                case 3:
                    surface[q].N[i][j]=F_N4(a,b);break;
             }
        }
}

void ElementUniversal::oblicanie_tabN() {
    for (int q = 0; q < 4; ++q) {
        int z = sqrt(E)-1;
        switch (q) {
            case 0:
                for (int i = 0; i < sqrt(E); ++i) {
                    funkcja_do_tabN(i,obj.x[i],-1,q);
                }break;
            case 1:
                for (int i = 0; i < sqrt(E); ++i) {
                    funkcja_do_tabN(i,1,obj.x[i],q);
                }break;
            case 2:
                for (int i = 0; i < sqrt(E); ++i) {
                    funkcja_do_tabN(i,obj.x[z],1,q); z--;
                }break;
            case 3:
                for (int i = 0; i < sqrt(E); ++i) {
                    funkcja_do_tabN(i,-1,obj.x[z],q); z--;
                }break;
        }
    }
}

void ElementUniversal::macierz_Hbc(int e,Element& elem,int nr_e) {//liczy tez wektor P
    double hbc_sciany[4][4];    fill_n(&hbc_sciany[0][0], 4 * 4, 0);
    double hbc_suma[4][4];      fill_n(&hbc_suma[0][0], 4 * 4, 0);
    double tmp[4][4];           fill_n(&tmp[0][0], 4 * 4, 0);
    double tmp_p[4];            fill_n(&tmp_p[0], 4, 0);
    double tmp_p_suma[4];       fill_n(&tmp_p_suma[0], 4, 0);
    double tmp_p_sciany[4];     fill_n(&tmp_p_sciany[0], 4, 0);

    for (int i = 0; i < 4; ++i) {
        double detj=0.;
        int z = (i == 3) ? 0 : i + 1;/*int z=i+1;  if(i==3){z=0;}*/
        if (elem.ID[i].BC && elem.ID[z].BC){
            detj =sqrt( pow((elem.ID[z].X - elem.ID[i].X),2) + pow((elem.ID[z].Y - elem.ID[i].Y),2) ) /2;
        }
                for (int l = 0; l < e; ++l) {
                    for (int j = 0; j < 4; ++j) {
                        for (int k = 0; k < 4; ++k) {
                            tmp[j][k] = surface[i].N[l][j] * surface[i].N[l][k];
                            tmp[j][k]*=obj.w[l];
                            hbc_suma[j][k]+=tmp[j][k];
                        }
                        tmp_p[j]=surface[i].N[l][j]* obj.w[l];
                        tmp_p_suma[j]+=tmp_p[j];
                    }
                }
                for (int j = 0; j < 4; ++j) {
                    for (int k = 0; k < 4; ++k) {
                        hbc_suma[j][k]*=(detj*dane.Alfa);
                        hbc_sciany[j][k]+=hbc_suma[j][k];
                        hbc_suma[j][k]=0;
                    }
                    tmp_p_suma[j]*=(dane.Tot*detj*dane.Alfa);
                    tmp_p_sciany[j]+=tmp_p_suma[j];
                    tmp_p_suma[j]=0;
                }

    }

    cout<<"HBC dla elementu "<<nr_e+1<<endl;
    display_tab_44(hbc_sciany);
    elem.set_Hbc(hbc_sciany);
    cout<<"\nWektor P dla elementu "<<nr_e+1<<endl;
    display_tab_4(tmp_p_sciany);
    elem.set_P(tmp_p_sciany);

}

void ElementUniversal::display_tab_22(double tab[2][2]){
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            cout<<setw(10)<<tab[i][j]<<"     ";
        }
        cout<<endl;
    }
}

void ElementUniversal::display_tab_4(double tab[4]){
    for (int i = 0; i < 4; ++i) {
            cout<<setw(10)<<setprecision(7)<<tab[i]<<"     ";
        }
        cout<<endl;
}

void ElementUniversal::display_tab_44(double tab[4][4]){
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            cout<<setw(10)<<tab[i][j]<<"     ";
        }
        cout<<endl;
    }
}
void ElementUniversal::display_tab_E4(double **tab){
    for (int i = 0; i < E; ++i) {
        for (int j = 0; j < 4; ++j) {
            cout<<setw(10)<<tab[i][j]<<"     ";
        }
        cout<<endl;
    }
    cout<<endl;
}
void ElementUniversal::display_H(double tab[4][4]){
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            cout << setw(10)<<tab[i][j] << "   ";
        }
        cout << endl;
    }
}

ElementUniversal::~ElementUniversal(){
    for (int i = 0; i < E; i++) {
        delete[] dtksi[i];
        delete[] dtetta[i];
    }
    delete[] dtksi;
    delete[] dtetta;
}

/*
Surface::~Surface() {
    for (int i = 0; i < e; ++i) {
        delete[] N[i];
    }
    delete[] N;
}*/
