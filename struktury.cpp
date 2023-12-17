#include "struktury.h"
#include <fstream>

void zapis_ParaView(ofstream& file, double *tab, Grid& siatka);

void Element::set_tabH(double tab[4][4]) {
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            tab_H[i][j] = tab[i][j];
        }
    }
}

void Element::set_Hbc(double tab[4][4]) {
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            Hbc[i][j] = tab[i][j];
        }
    }
}
void Element::set_C(double tab[4][4]) {
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            C[i][j] = tab[i][j];
        }
    }
}
void Element::set_P(double tab[4]) {
    for (int i = 0; i < 4; ++i) {
        P[i] = tab[i];
    }
}

void Element::obliczanie_h_calk() {
    // double H_calk[4][4];  fill_n(&H_calk[0][0], 4 * 4, 0);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            H_CALK[i][j] = tab_H[i][j] + Hbc[i][j];
        }
    }
}
void Element::display_H_calk(){
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            cout << setw(10)<<H_CALK[i][j] << "   ";
        }
        cout << endl;
    }
}
void Element::display_P_calk(){
    for (int i = 0; i < 4; ++i) {
            cout << setw(10)<<P[i] << "   ";
    } cout << endl;
}

void Element::display_C_calk(){
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            cout << setw(10)<<C[i][j] << "   ";
        }
        cout << endl;
    }
}

Grid::Grid(){
    Nodes = new Node[dane.Nodes_number];
    Elements = new Element[dane.Elements_number];
    for (int i = 0; i < dane.Nodes_number; i++) {Nodes[i] = Node();}
    //      for (int x = 0; x < dane.Elements_number; x++) {Elements[x] = Node();}}
    // ~Grid() { delete[] Nodes; delete[] Elements; }
};

SOE::SOE(){

   G_H = new double*[dane.Nodes_number];
   G_P = new double [dane.Nodes_number];
   G_C = new double*[dane.Nodes_number];

    C_tau = new double*[dane.Nodes_number];
    t0 = new double [dane.Nodes_number];
    C_t0 = new double [dane.Nodes_number];

    for (int i = 0; i < dane.Nodes_number; ++i) {
        G_H[i] = new double[dane.Nodes_number];
        G_C[i] = new double[dane.Nodes_number];
        C_tau[i] = new double[dane.Nodes_number];
    }

    for (int i = 0; i < dane.Nodes_number; ++i) {
        G_P[i]=0.0;
        C_t0[i]=0.0;
        t0[i]=dane.InitialTemp;
        for (int j = 0; j < dane.Nodes_number; ++j) {
            G_H[i][j]=0.0;
            G_C[i][j]=0.0;
            C_tau[i][j]=0.0;
        }
    }

};

void SOE::agregacja(Grid &siatka) {
    for (int i = 0; i < dane.Elements_number; ++i) {
        for (int k = 0; k < 4; ++k) {
            for (int l = 0; l < 4; ++l) {
                int a,b;
                    a = siatka.Elements[i].ID[k].id-1;
                    b = siatka.Elements[i].ID[l].id-1;
                G_H[a][b]+=siatka.Elements[i].H_CALK[k][l];
                G_C[a][b]+=siatka.Elements[i].C[k][l];
            }
        }

        for (int j = 0; j < 4; ++j) {
            int c;
            c=siatka.Elements[i].ID[j].id-1;
            G_P[c]+=siatka.Elements[i].P[j];
        }
    }

}

//uklad rownan - metoda 1
/*
// wyswietl macierz
void SOE::print(double ** A)
{
    for (int i=0; i<dane.Nodes_number; i++, cout<<endl)
        for (int j=0; j<dane.Nodes_number; j++)
            cout<<" | "<<setw(10)<<A[i][j];

    cout<<endl;
}*/

// zamien dwa rzedy ze soba
void SOE::swap_row(double **A, double * B, int i, int j)
{
    //cout<<"Zamieniam rzad "<< i <<" z "<<j<<endl;
    double temp2 = B[i];
    B[i] = B[j];
    B[j] = temp2;

    for (int k=0; k<dane.Nodes_number; k++)
    {
        double temp = A[i][k];
        A[i][k] = A[j][k];
        A[j][k] = temp;
    }
}

// redukuje macierz do postaci schodkowej
int SOE::eliminacja(double **A,  double *B)
{
    for (int k=0; k<dane.Nodes_number; k++)
    {
        // wartosc maksymalna i indeks piwotu
        int i_max = k;                      //indeks
        int v_max = A[i_max][k];          //wartosc max

        //znalezc maksymalna amplitude dla piwotu
        for (int i = k+1; i < dane.Nodes_number; i++)
            if (abs(A[i][k]) > v_max)
                v_max = A[i][k], i_max = i;

        // sprawdzenie warunku - jesli na przekatnej jest zero
        // to doprowadzi to pozniej do dzielenia przez zero
        if (!A[k][i_max])
            return k; // macierz jednostkowa

        // Zamiana rzedu maksymalnego z obecnym
        if (i_max != k)
            SOE::swap_row(A, B, k, i_max);

        for (int i=k+1; i<dane.Nodes_number; i++)
        {
            //f - kty element bierzacego rzedu na zero
            double f = A[i][k]/A[k][k];

            for (int j=k+1; j<=dane.Nodes_number; j++){
                if (j == dane.Nodes_number)
                    B[i] -= B[k]*f;
                else
                    A[i][j] -= A[k][j]*f; // odjecie f*odpowiedznie elemeny wiersza
            }
            //dolna czesc kolumny zapelniamy zerami
            A[i][k] = 0;
        }
        //print(A);        //etap obliczen
    }
    //SOE::print(A);            //macierz schodkowa
    return -1;
}

//obliczenie niewiadomych
double* SOE::result(double **A,  double *B, double *x, int q, Grid& siatka){

    // od ostatniego rownania do pierwszego
    for (int i = dane.Nodes_number-1; i >= 0; i--)
    {
        x[i] = B[i];// wspolczynniki macierzy B

        //j=i+1, dlatego ze macierz jest trojkatna gorna
        for (int j=i+1; j<dane.Nodes_number; j++)
        {
            x[i] -= A[i][j]*x[j]; //odejmujemy wszystkie wartosci oprocz wsp zmiennej aktualnie liczonej
        }
        //dzielimy wartosci po = przez wspolczynnik prz x i otrzymujemy szukany wynik
        x[i] = x[i]/A[i][i];
    }

    //zapis do pliku ParaView
    const string name = "foo";
    const string extension = ".vtk";

    string full_name = name + to_string(q) + extension;
    ofstream file_name(full_name);
    zapis_ParaView(file_name,x,siatka);
    //


//    cout<<endl<<endl<<"Rozwiazanie ukladu rownan:"<<endl;                                                          //1
    for (int i=0; i<dane.Nodes_number; i++){
        t0[i]=x[i];
//        cout<<setprecision(7)<< " x"<<i<<" = "<< x[i] <<endl;                                                      //1
    }
    double max=x[0],min=x[0];
    for (int i = 0; i < dane.Nodes_number; ++i) {
        if(x[i]>max){
            max =x[i];
        }
        if(x[i]<min){
            min=x[i];
        }
    }
    cout<<endl<<"Min: "<<min<<" Max: "<<max<<endl;                                                                 //1
    return x;
}

void SOE::roz_gauss(double *x,double **h,double *tp,int q,Grid& siatka)
{
    // redukcja do postaci schodkowej
    int singularM = SOE::eliminacja(h, tp);

    // w przpadku macierzy jednostkowej
    if (singularM != -1)
    {
        printf("Macierz jednostkowa!\n");

        if (tp[singularM])
            printf("Rownanie sprzeczne!");
        else
            printf("Nieskonczona liczba rozwiazan!");

        return;
    }

    result(h, tp, x,q,siatka);
}

//uklad rownan koniec

void SOE::rozwiazywanie_temp(Grid& siatka){
    double* X = new double[dane.Nodes_number];
    double dt = dane.SimulationStepTime;
    double time=500;
    double *tmp_G_P = new double[dane.Nodes_number];

    for (int i = 0; i < dane.Nodes_number; ++i) {
        for (int j = 0; j < dane.Nodes_number; ++j) {
            C_tau[i][j]=G_C[i][j]/dane.SimulationStepTime;    //([C]/tau)
            G_H[i][j]+=C_tau[i][j];                          // [H]+[C_tau] to jest [H]+([C]/tau)
        }
    }
    for (int i = 0; i < (time/dt); ++i) {
        for (int k = 0; k < dane.Nodes_number; ++k) {
            tmp_G_P[k]=G_P[k];
            C_t0[k]=0.0;
        }
        //kopia
        double **tmp_G_H = new double*[dane.Nodes_number];
        for (int w = 0; w < dane.Nodes_number; ++w) {
            tmp_G_H[w] = new double[dane.Nodes_number];
        }
        for (int r = 0; r < dane.Nodes_number; ++r) {
            for (int j = 0; j < dane.Nodes_number; ++j) {
                tmp_G_H[r][j] = G_H[r][j];
            }
        }//kopia
        for (int q = 0; q < dane.Nodes_number; ++q) {
            for (int j = 0; j < dane.Nodes_number; ++j) {
                C_t0[q]+= (C_tau[q][j] * t0[j]);                    //([C]/tau)*t0   //t0 bedzie sie zmieniac - kolejne rozwiazania temp
            }
             tmp_G_P[q]+=C_t0[q]; // [C_t0]+[P]          //   ([C]/tau)*t0 + [P]
        }

        roz_gauss(X,tmp_G_H,tmp_G_P,i,siatka);
    }

}



void SOE::display_G_H(){
    for (int i = 0; i < dane.Nodes_number; ++i) {
        for (int j = 0; j < dane.Nodes_number; ++j) {
            cout <<setprecision(5)<<G_H[i][j] << "   ";
        }
        cout << endl;
    }
}
void SOE::display_G_P(){
    for (int i = 0; i < dane.Nodes_number; ++i) {
            cout <<setprecision(6)<<G_P[i] << "   ";
        }
    cout << endl;
}
void SOE::display_G_C(){
    for (int i = 0; i < dane.Nodes_number; ++i) {
        for (int j = 0; j < dane.Nodes_number; ++j) {
            cout <<setprecision(5)<<G_C[i][j] << "   ";
        }
        cout << endl;
    }
    cout << endl;
}

/*SOE::~SOE() {
   for (int i = 0; i < dane.Nodes_number; ++i) {
        delete[] G_H[i];
    }
    delete[] G_H;
}*/

void wyswietl_dane(){
    cout<< dane.SimulationTime <<endl<< dane.SimulationStepTime<<endl
        << dane.Conductivity<<endl<< dane.Alfa<<endl<< dane.Tot<<endl
        << dane.InitialTemp<<endl<< dane.Density<<endl<< dane.SpecificHeat<<endl
        << dane.Nodes_number<<endl<< dane.Elements_number<<endl;
}

void wyswietl_nodes(Grid siatka){
    for(int i=0;i<dane.Nodes_number;i++){
        cout<<"id: "<<setw(3)<<siatka.Nodes[i].id<<"      X[ "
            <<setw(13)<<setprecision(12)<<siatka.Nodes[i].X<<" ]       Y[ "
            <<setw(13)<<setprecision(12)<<siatka.Nodes[i].Y<<" ]"<<" Status "<<siatka.Nodes[i].BC<<endl;
    }
}
void wyswietl_elements(Grid siatka){
    for (int i = 0; i <dane.Elements_number; ++i) {
        cout<<endl<<"Element "<<i+1<<": "<<endl;
        for (auto & j : siatka.Elements[i].ID) { // for (int j = 0; j < 4; ++j) {
            cout<<"["<<j.X<<" , ";                 //       cout<<siatka.Elements[i].ID[j].X<<" ";
            cout<<j.Y<<"] BC: ";
        cout<<j.BC<<endl;}                       //       cout<<siatka.Elements[i].ID[j].Y<<" | ";}
        cout<<endl;
    }
}

void wyswietl_BC(Grid siatka){
    for (int j = 0; j < dane.Nodes_number; ++j) {
        cout<<setw(3)<<siatka.Nodes[j].id<<" ";
    }cout<<endl;
    for (int j = 0; j < dane.Nodes_number; ++j) {
        cout<<setw(3)<<siatka.Nodes[j].BC<<" ";
    }
}
