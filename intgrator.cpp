//============================================================================
// Name        : Integration
// Author      : C. Poisson
// Version     : 1.1
// Description :
//============================================================================

#include "integrator.h"
#include <cmath>

//============================================================================
// Name: integrate_trapz
// Parameter f: A mathematical function
// Parameter a: Lower bound
// Parameter b: Upper bound
// Description: Compute the integral of the function  f in [a,b] using the
// trapeze method
// Return I: The value of the integral
//============================================================================

double integrate_trapz(double(*f)(double), double a, double b){
	double J=1;
	double I=(b-a)*(f(a)/2.+f(b)/2.);// First itteration with 1 trapeze
	double oldI=0.;
	double eps=1E-8; // Setting the precision

	//Trapezoidal rule
	//Integral f = Δx/2 *(f(x0)+2f(x1)+...+2f(xn-1)+f(xn))
	do{
		oldI = I;// Stocking old value
		J++; // Increamenting
		double h=(b-a)/pow(2.,J-1.); // New interval length
		double s=0.;
		double x=a+h;// The new first abscissa

		for(double k=0; k<=pow(2.,J-1)-2.; k+=2.){
			s+=f(x);
			x+=2.*h;
		}

		I=0.5*I+h*s; // Computing the integral
			//		I1=(b-a)*(f(a)/2.+f(b)/2.)
			//		I2=(b-a)/2 *(f(a)/2.+f(c)/2.)+(b-a)/2 *(f(b)/2. + f(c)/2.)
			// <=>	I2=I1/2+(b-a)/2 *f(c)
			// <=>	I2=I1/2+h/2 *f(c)
			// <=>	I2=0.5*(I+h*f(c))
	} while(fabs(oldI-I)>eps);// Checking precision

    if(I<1E-10){
        return 0;
    }
    else{
        return I;
    }

}

//============================================================================
// Name: integrate_trapz
// Parameter f: A mathematical function
// Parameter a: Lower bound
// Parameter b: Upper bound
// Parameter p: Function parameter
// Description: Compute the integral of the parametric function f(x,p) in [a,b]
// using the trapeze method
// Return I: The value of the integral
//============================================================================

double integrate_trapz(double(*f)(double, double), double a, double b, double p){
	double J=1;
	double I=(b-a)*(f(a,p)/2.+f(b,p)/2.);// First itteration with 1 trapeze
	double oldI=0.;
	double eps=1E-6; // Setting the precision

	//Integral f = Δx/2 *(f(x0)+2f(x1)+...+2f(xn-1)+f(xn))
	do{
		oldI=I;// Stocking old value
		J++; // Increamenting
		double h=(b-a)/pow(2.,J-1.); // New interval length
		double s=0.;
		double x=a+h;// The new first abscissa

		for (double k=0; k<=pow(2.,J-1)-2.; k+=2.){
			s+=f(x,p);
			x+=2.*h;
		}

		I=0.5*I+h*s; // Computing the integral
			//		I1=(b-a)*(f(a)/2.+f(b)/2.)
			//		I2=(b-a)/2 *(f(a)/2.+f(c)/2.)+(b-a)/2 *(f(b)/2. + f(c)/2.)
			// <=>	I2=I1/2+(b-a)/2 *f(c)
			// <=>	I2=I1/2+h/2 *f(c)
			// <=>	I2=0.5*(I+h*f(c))
	} while (fabs(oldI-I)>eps);// Checking precision

	if(I<1E-10){
        return 0;
    }
    else{
        return I;
    }
}

//============================================================================
// Name: integrate
// Parameter f: A mathematical function
// Parameter a: Lower bound
// Parameter b: Upper bound
// Description: Compute the integral of the parametric function f(x) in [a,b]
// using the simpson method
// Notes: Converges faster than trapeze method
// Return I: The value of the integral
//============================================================================

double integrate(double(*f)(double), double a, double b){
	int J=2;
	double oldS =0.;
	double oldI=(b-a)*(f(a)/2.+f(b)/2); // First itteration
	double I=0.5*oldI+(b-a)/2.*f((a+b)/2.);// Second itteration
	double S = 4./3.*I-1./3.*oldI;
	//
	double eps=1E-6;

	// Simpson's rule
	// Integral f = Δx/3 *(f(x0)+4f(x1)+2f(x2)...+2f(xn-2)+4f(xn-1)+f(xn))
	do{
		// Stocking old values
		oldS=S;
		oldI=I;
		J++;// Increamenting
		double h=(b-a)/pow(2.,J-1); // New interval length
		double s=0.;
		double x=a+h; // The new first abscissa

		for (double k=0; k<=pow(2.,J-1)-2.; k+=2.){
			s=s+f(x);
			x=x+2.*h;
		}

		I=I*0.5+h*s;
		S=4./3.*I-1./3.*oldI;


	} while (fabs(oldS-S)>eps);

    if(S<1E-10){
        return 0;
    }
    else{
        return S;
    }
}

double integrate(double(*f)(double, double), double a, double b, double p){
	int J=2;
	double oldS=0.;
	double oldI=(b-a)*(f(a,p)/2.+f(b,p)/2.);
	double I = 0.5*oldI+(b-a)/2. *f((a+b)/2.,p);
	double S = 4./3. *I-1./3.*oldI;
	double eps = 0.000000001;

	do{
		oldS=S;
		oldI=I;
		J++;
		double h=(b-a)/pow(2.,J-1);
		double s=0.;
		double x=a+h;

		for (double k=0;k<=pow(2.,J-1)-2.; k+=2.){
			s=s+f(x,p);
			x=x+2.*h;
		}

		I = I*0.5+h*s;
		S = 4./ 3.*I -1./3.*oldI;

	} while (fabs(oldS-S)>eps);

	if(S<1E-10){
        return 0;
    }
    else{
        return S;
    }
}

