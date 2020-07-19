//
//  Question3.cpp
//  MA323 Project
//
//  Created by Tyler Mitchell on 03/06/2020.
//  Copyright © 2020 Tyler Mitchell. All rights reserved.
//



// Question3.cpp
#include "Question3.hpp"
#include "Normals.hpp"
#include "RandomNumberGenerators.hpp"

// QUESTION 3 CODE //

double analytical(double C, double H, double S0, double r, double T, double sigma){

    double d = ( log((S0/H)) + (r-(0.5*sigma*sigma))*T )/(sigma*sqrt(T));
    double result = C*exp(-1*r*T)*CumulativeNormal(d);
    return result; // The analytical price of our digital option.
}

double Question3Base::S_T(double S0,double r,double T,double sigma){

            double Z, diffusion, drift = (r-(0.5*sigma*sigma))*T;
            Z = normal_rv(0.0,1.0);
            diffusion = sigma*sqrt(T)*Z;
            return (S0*exp(drift + diffusion)); // Our simulated GBM stock price
}

MonteCarlo::MonteCarlo(double C, double H, double S0, double r, double T, double sigma, int n){
    // Initialising private member variables.
        myC     = C;
        myH     = H;
        myS0    = S0;
        myr     = r;
        myT     = T;
        mysigma = sigma;
        myn     = n;

}

double MonteCarlo::getEstimate(void){

    double ST, sum = 0;
    for(int i = 0; i < myn; i++){
        // Performing myn = n Monte-Carlo iterations, to produce a numerical estimate of the
        // digital option price.

        // Generating Our stock path.
        ST = S_T(myS0,myr,myT,mysigma);

        // Our Indicator Function.
        if( ST > myH ){
            sum += (myC*exp(-1*myr*myT));
        }
    }

    myEstimate = sum/myn; // Our Monte-Carlo Estimate.


    // Calculating Option analytical Variance.
    double d = ( log(myS0/myH) + (myr-(0.5*mysigma*mysigma))*myT )/(mysigma*sqrt(myT));
    myanalyticalVariance = myC*myC*exp(-2*myr*myT)*CumulativeNormal(d)*CumulativeNormal(-1*d)/myn; // assigning value to respective private member variable.

    return myEstimate;
}

double MonteCarlo::ConfidenceInterval(double epsilon){
    // Calculating Confidence Interval.
    // Function returns the interval width.

    double alpha = InverseCumulativeNormal(1.0 - (epsilon/2.0));
    myleftCI = myEstimate - alpha*sqrt(myanalyticalVariance);  // assigning values to respective
    myrightCI = myEstimate + alpha*sqrt(myanalyticalVariance); // private member variables.
    return myrightCI - myleftCI; // The CI interval width.
}

ControlVariates::ControlVariates(double C, double H, double S0, double r, double T, double sigma, int n){
    // Initialising private member variables.
        myC = C;
        myH = H;
        myS0 = S0;
        myr = r;
        myT = T;
        mysigma = sigma;
        myn = n;
}

double ControlVariates::bstarhat(void){

    double numeratorSum = 0, denominatorSum = 0;
    double Xi,Yi;

    // Calculating EX and EY =(approx)= Monte-Carlo(Y).
    double EX = exp(myr*myT)*myS0;
    MonteCarlo MC(myC,myH,myS0,myr,myT,mysigma,myn);
    double EY = MC.getEstimate();

    for(int i = 0; i < myn; i++){
        Xi = S_T(myS0,myr,myT,mysigma);
        Yi = 0;
        if( Xi > myH ){
            Yi = myC*exp(-1*myr*myT);
        }
        numeratorSum += ((Xi - EX)*(Yi - EY));
        denominatorSum += ((Xi - EX)*(Xi - EX));
    }
    return numeratorSum/denominatorSum; // Our B*hat estimate.

}

double ControlVariates::getEstimate(void){

    double yib, ST, CVEstimate, sum = 0, squaresum = 0;
    double yi,Xi, EX = exp(myr*myT)*myS0; // E[X] = E[ST].
    double b = bstarhat();                // Our b* estimate.
    for(int i = 0; i < myn; i++){

        // Generating Our stock path.
        ST = S_T(myS0,myr,myT,mysigma);

        // Calculating Yi.
        yi = 0;
        if( ST > myH ){
            yi = myC*exp(-1*myr*myT);
        }

        // Calculating Yi(b).
        Xi = ST;
        yib = yi - b*(Xi - EX);

        sum += yib;
        squaresum += (yib*yib);

    }

    // Calculating statistics.
    CVEstimate = sum/myn;
    myvariance = (squaresum - (myn*CVEstimate*CVEstimate))/(myn*(myn-1)); // Our sample variance estimate of our Control Variates estimator.

    return CVEstimate; // Our Control Variates Estimate.
}
