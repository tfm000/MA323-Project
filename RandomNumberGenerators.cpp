//
//  RandomNumberGenerators.cpp
//  MA323 Project
//
//  Created by Tyler Mitchell on 03/06/2020.
//  Copyright © 2020 Tyler Mitchell. All rights reserved.
//

#include "RandomNumberGenerators.hpp"
#include<cstdlib>
#include<cmath>

double myuniform(void){
    int myrand = rand();
    while((myrand==0)||(myrand==RAND_MAX)){myrand = rand();}
    double myuni = myrand/static_cast<double>(RAND_MAX);
    return myuni;
}

double uniformab(double a,double b){
    double u = myuniform();
    return a + ((b-a)*u);
}

double exponential_rv(double lambda){
    double u = myuniform();
    return -1*log(u)/lambda;
}

double normal_rv(double mu, double sigma){
    double myuni1 = myuniform(), myuni2 = myuniform();
    double myexp = -2.0*log(myuni1);
    double mytheta = 2.0*M_PI*myuni2;
    double normal = mu + sigma*(sqrt(myexp)*cos(mytheta));
    return normal;
}


