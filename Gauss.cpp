#include "Gauss.h"

Gauss::Gauss (int N): n(N){  //punktowanie ile n
    x=new double[n]; //pkt - x
    w=new double[n]; // wagi
    switch (N) {
        case 1: x[0]=0.0; w[0]=2.; break;
        case 2: x[0] = -sqrt(1./3.); x[1] = sqrt(1./3.);  w[0] = w[1] = 1.0;break;
        case 3: x[0] = -sqrt(3./5.); x[1] = 0.0; x[2] = sqrt(3./5.);
            w[0] = w[2] = 5.0 / 9.0; w[1] = 8.0 / 9.0;break;
        case 4:
            x[0] = -sqrt( 3./7. + 2./7.*sqrt(6./5.));
            x[1] = -sqrt( 3./7. - 2./7.*sqrt(6./5.));
            x[2] = sqrt( 3./7. - 2./7.*sqrt(6./5.));
            x[3] = sqrt( 3./7. + 2./7.*sqrt(6./5.));
            w[0] = w[3] = (18. - sqrt(30.)) / 36.;
            w[1] = w[2] = (18. + sqrt(30.)) / 36.; break;
        case 5:
            x[0] = 1./3. * sqrt(5. + 2. * sqrt(10./7.));
            x[1] = - 1./3. * sqrt(5. - 2. * sqrt(10./7.));
            x[2] = 0.0;
            x[3] = 1./3. * sqrt(5. - 2. * sqrt(10./7.));
            x[4] = - 1./3. * sqrt(5. + 2. * sqrt(10./7.));
            w[0] = w[4] = (322. - 13.*sqrt(70.)) / 900.;
            w[1] = w[3] = (322. + 13.*sqrt(70.)) / 900.;
            w[2] = 128. / 225.; break;
        default: exit;
    }

}

void Gauss::wym1D(double(*f)(double)){
    double result;
    for(int i=0;i<n;i++){result+=w[i]*f(x[i]);}
    cout<<result<<endl;
}

void Gauss::wym2D(double(*f)(double,double)){
    double result;
    for(int i=0;i<n;i++){
        for (int j = 0; j < n; j++) {
            result+=w[j]*w[i]*f(x[i], x[j]);
        }}
    cout<<result<<endl;

}

Gauss::~Gauss(){
    delete[] w;
    delete[] x;
}

