#include <iostream>
#include <cmath>
#include "integrator.h"



using namespace std;

double f(double x){
    return x*x;
}
double f2(double x){
    return cos(x);
}
double f3(double x){
    return exp(x);
}


int main(){
    cout.precision(6);
    cout<< "Function:x*x [0;3]\t\t"<< "Trapezoidal method: "<<integrate_trapz(f,0,3)<<"\t\tSimpson's method: "<< integrate(f,0,3)<< endl;
    cout<< "Function:Cos(x)[0;Pi]\t\t"<< "Trapezoidal method: "<<integrate_trapz(f2,0,M_PI)<<"\t\tSimpson's method: "<< integrate(f2,0,M_PI)<< endl;
    cout<< "Function:Cos(x)[0;Pi/2]\t\t"<< "Trapezoidal method: "<<integrate_trapz(f2,0,M_PI_2)<<"\t\tSimpson's method: "<< integrate(f2,0,M_PI_2)<< endl;
    cout<< "Function:Exp(x)[0;1]\t\t"<< "Trapezoidal method: "<<integrate_trapz(f3,0,1)<<"\tSimpson's method: "<< integrate(f3,0,1)<< endl;


}







