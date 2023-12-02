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
    cout<<"DANE ile node: "<<dane.Nodes_number<<endl;
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

    display_G_H();
    cout<<"P"<<endl;
    display_G_P();cout<<endl;
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
