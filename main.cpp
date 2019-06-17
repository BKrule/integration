#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <fstream>
#include "integrator.h"



using namespace std;

double f(double x) {
	return x * x;
}
double f2(double x) {
	return cos(x);
}
double f3(double x) {
	return exp(x);
}
double S(double t) {
	return sin(t*t);
}
double C(double t) {
	return cos(t*t);
}



int main() {
	cout.precision(8);
	cout << "Function:x*x [0;3]\t\t" << "Trapezoidal method: " << integrate_trapz(f, 0, 3) << "\t\tSimpson's method: " << integrate_simps(f, 0,3) << endl;
	cout << "Function:Cos(x)[0;Pi]\t\t" << "Trapezoidal method: " << integrate_trapz(f2, 0, M_PI) << "\t\tSimpson's method: " << integrate_simps(f2, 0, M_PI) << endl;
	cout << "Function:Cos(x)[0;Pi/2]\t\t" << "Trapezoidal method: " << integrate_trapz(f2, 0, M_PI_2) << "\t\tSimpson's method: " << integrate_simps(f2, 0, M_PI_2) << endl;
	cout << "Function:Exp(x)[0;1]\t\t" << "Trapezoidal method: " << integrate_trapz(f3, 0, 1) << "\tSimpson's method: " << integrate_simps(f3, 0, 1) << endl;

	ofstream outS("S(x).data", ios::trunc);
	if (outS.is_open()) {
		for (double c=0; c<10; c+=0.01) {
			outS<<c<<"\t"<<integrate_simps(S, 0, c)<<endl;
		
			
		}
	}
	else {
		cout<<"Output error"<<endl;
	}
	ofstream outC("C(x).data", ios::trunc);
	if (outC.is_open()) {
		for (double c=0; c<10; c+=0.01) {
			outC<<c<<"\t"<<integrate_simps(C, 0, c)<<endl;


		}
	}
	else {
		cout<<"Output error"<<endl;
	}

	return 0;
}







