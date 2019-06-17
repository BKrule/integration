#ifndef integrator_hpp
#define integrator_hpp

#include <stdio.h>

double integrate_trapz(double (*f)(double), double a, double b);

double integrate_trapz(double (*f)(double, double), double a, double b, double p);

double integrate_simps(double (*f)(double), double a, double b);

double integrate_simps(double (*f)(double, double), double a, double b, double p);



#endif /* integrator_hpp */
#pragma once
