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
Grid::Grid(){
    Nodes = new Node[dane.Nodes_number];
    Elements = new Element[dane.Elements_number];
    for (int i = 0; i < dane.Nodes_number; i++) {Nodes[i] = Node();}
    //      for (int x = 0; x < dane.Elements_number; x++) {Elements[x] = Node();}}
    // ~Grid() { delete[] Nodes; delete[] Elements; }
};

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
