#include "struktury.h"

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
void Element::set_P(double tab[4]) {
    for(int i = 0; i < 4; i++) {
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
    for (int i = 0; i < dane.Nodes_number; ++i) {
        G_H[i] = new double[dane.Nodes_number];
    }
    for (int i = 0; i < dane.Nodes_number; ++i) {
        G_P[i]=0.0;
        for (int j = 0; j < dane.Nodes_number; ++j) {
            G_H[i][j]=0.0;
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
            }
        }

        for (int j = 0; j < 4; ++j) {
            int c;
            c=siatka.Elements[i].ID[j].id-1;
            G_P[c]+=siatka.Elements[i].P[j];
        }
    }
}

//uklad rownan
//metoda druga
// wyswietl macierz
void SOE::print(double ** A)
{
    for (int i=0; i<dane.Nodes_number; i++, cout<<endl)
        for (int j=0; j<dane.Nodes_number; j++)
            cout<<" | "<<setw(10)<<A[i][j];

    cout<<endl;
}

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
double* SOE::result(double **A,  double *B, double *x){

    // od istatniego rownania do pierwszego
    for (int i = dane.Nodes_number-1; i >= 0; i--)
    {
        // wspolczynniki macierzy B
        x[i] = B[i];

        //j=i+1, dlatego ze macierz jest trojkatna gorna
        for (int j=i+1; j<dane.Nodes_number; j++)
        {
            //odejmujemy wszystkie wartosci oprocz wsp zmiennej aktualnie liczonej
            x[i] -= A[i][j]*x[j];
        }

        //dzielimy wartosci po = przez wspolczynnik prz x i otrzymujemy szukany wynik
        x[i] = x[i]/A[i][i];
    }

    cout<<endl<<"Rozwiazanie ukladu rownan:"<<endl;
    for (int i=0; i<dane.Nodes_number; i++)
        cout<<setprecision(7)<< " x("<<i<<") = "<< x[i] <<endl;
    return x;
}

void SOE::gauss_crout(double *x)
{
    // redukcja do postaci schodkowej
    int singularM = SOE::eliminacja(G_H, G_P);

    // w przpadku macierzy jednostkowej
    if (singularM != -1)
    {
        printf("Macierz jednostkowa!\n");

        if (G_P[singularM])
            printf("Rownanie sprzeczne!");
        else
            printf("Nieskonczona liczba rozwiazan!");

        return;
    }

    result(G_H, G_P, x);
    cout<<endl;
}

//metoda trzecia

bool SOE::gauss_crout2(double *X)
{
    double **AB = new double *[dane.Nodes_number];   //tablica dwuwymiarowa do obliczen
    for (int i = 0; i <= dane.Nodes_number; i++)
        AB[i] = new double[dane.Nodes_number+1];

    for(int i = 0; i < dane.Nodes_number ;i++)
        for(int k = 0; k <= dane.Nodes_number ;k++)
            if(k == dane.Nodes_number)
                AB[i][k] = G_P[i];
            else
                AB[i][k] = G_H[i][k];

    //pomocnicza tablica zamiany kolumn
    int *W = new int[dane.Nodes_number];
    for(int i=0; i<dane.Nodes_number ;i++)
        W[i] = i;

    //iteratory petli
    int    i, j, k;
    double m, s;        //zmienne pomocnicze

    for( i = 0; i <= dane.Nodes_number; i++ )
        W [i] = i;

    // eliminacja współczynników
    for( i = 0; i < dane.Nodes_number - 1; i++ )
    {
        k = i;
        for( j = i + 1; j < dane.Nodes_number; j++ )
            if( fabs ( AB[i][ W[k] ] ) < fabs ( AB[i][ W[j] ] ))
                k = j;
        swap ( W [k], W [i] );


        for( j = i + 1; j < dane.Nodes_number ;j++ )
        {
            if( fabs (AB[i][ W[i] ]) < eps ) return false;
            m = -AB[j][ W [i] ] / AB[i][ W[i] ];
            for( k = i + 1; k <= dane.Nodes_number; k++ )
                AB[j][ W[k] ] += m * AB[i][ W[k] ];
        }
    }

    // wyliczanie niewiadomych
    for( i = dane.Nodes_number - 1; i >= 0; i--)
    {
        if( fabs (AB[i][ W[i] ]) < eps ) return false;
        s = AB[i][dane.Nodes_number];
        for( j = dane.Nodes_number - 1; j >= i + 1; j--)
            s -= AB[i][ W[j] ] * X [ W[j] ];
        X [ W[i] ] = s / AB[i][ W[i] ];
    }

    //print(AB);

    for( i = 0; i < dane.Nodes_number; i++ )
        cout << "x" << i + 1 << " = " << setw ( 7 ) << X [i]
             << endl;

    //for (int i = 0; i < N; i++)
    //delete [] AB[i];
    //delete [] AB;

    return true;
}

//uklad rownan

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
        cout<<"Element "<<i+1<<": "<<endl;
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
