#include <fstream>
#include <string>
#include <sstream>

#include "ElementUniversal.h"

GlobalData dane;

void wczytywanie(Grid& siatka){
    string zm,sm; char smiec; int liczba, liczba1, l1,l2,l3,l4; double w_x,w_y;
    fstream plik;
   // plik.open("Test1_4_4.txt", ios::in);
    plik.open("Test2_4_4_MixGrid.txt", ios::in);
    //plik.open("Test3_31_31_kwadrat.txt", ios::in);
    //plik.open("Testowe.txt", ios::in);

    if(plik.is_open())//plik.good())
    {
        for(int i=0;i<10;i++){
            try {
                //getline(plik, zm,' ');  plik>>liczba;
                getline(plik, zm);  istringstream jedna(zm);
                jedna>>sm>>liczba;
                switch(i){
                    case 0:dane.SimulationTime=liczba;break;
                    case 1:dane.SimulationStepTime=liczba;break;
                    case 2:dane.Conductivity=liczba;break;
                    case 3:dane.Alfa=liczba;break;
                    case 4:dane.Tot=liczba;break;
                    case 5:dane.InitialTemp=liczba;break;
                    case 6:dane.Density=liczba;break;
                    case 7:dane.SpecificHeat=liczba;break;
                    case 8:dane.Nodes_number=liczba;break;
                    case 9:dane.Elements_number=liczba;break;
                    default:break;
                }
            } catch (const invalid_argument& e) {
                cerr << "Blad konwersji na liczbe calkowita " << e.what() << endl;
            }
        }
        cout<<"dane:"<<endl;wyswietl_dane();
        //Grid siatka;
        siatka=Grid();
        getline(plik, zm);//"nieuzyteczna linijka z pliku;
        for(int i=0;i<dane.Nodes_number;i++){
            try {
                getline(plik, zm);  istringstream linia(zm);
                linia>>liczba>>smiec>>w_x>>smiec>>w_y;
                siatka.Nodes[i].id=liczba;
                siatka.Nodes[i].X=w_x;
                siatka.Nodes[i].Y=w_y;
            } catch (const invalid_argument& e) {
                cerr << "Blad " << e.what() << endl;
            }
        }
        cout<<"wspol node:"<<endl;wyswietl_nodes(siatka);
        getline(plik, zm);//nieuztyczena linika z pliku";
        for(int i=0;i<dane.Elements_number;i++){
            try {
                getline(plik, zm);  istringstream linia1(zm);
                linia1>>liczba1>>smiec>>l1>>smiec>>l2>>smiec>>l3>>smiec>>l4;
                siatka.Elements[i].ide=liczba1;
                siatka.Elements[i].ID[0].id=l1;
                siatka.Elements[i].ID[1].id=l2;
                siatka.Elements[i].ID[2].id=l3;
                siatka.Elements[i].ID[3].id=l4;

            } catch (const invalid_argument& e) {
                cerr << "Blad " << e.what() << endl;
            }
        }
       // cout<<"elem :"<<endl;//wyswietl_elements(siatka);
       getline(plik, zm);//cout<<"nieuzyteczna linijka z pliku :"<<zm<<endl;
        getline(plik, zm);  istringstream linijka(zm) ;
       int tmp;
        while(linijka>>tmp) {
            linijka.ignore();//igonruje koljeny znak - przecinki
             siatka.Nodes[tmp-1].BC=1;
        }
        cout<<endl<<"BC:"<<endl;wyswietl_BC(siatka);
        
        //przypisanie wszytskiego
        for (int i = 0; i < dane.Elements_number; ++i){
            for (int j = 0; j < 4; ++j) {
                siatka.Elements[i].ID[j]=siatka.Nodes[siatka.Elements[i].ID[j].id-1];
            }
        }

        cout<<endl;
        wyswietl_elements(siatka);

        /*for(auto & node : grid.elements[i].nodes){
            node = grid.nodes[node.node_id-1];
        }*/
            
        cout<<endl;
        plik.close();
    }else{
        cerr << "Nie udalo się otworzyc pliku." << endl;
    }

}
