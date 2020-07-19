//
//   Question4.cpp
//  MA323 Project
//
//  Created by Tyler Mitchell on 03/06/2020.
//  Copyright © 2020 Tyler Mitchell. All rights reserved.
//

#include "Question4.hpp"
#include "Normals.hpp"
#include "RandomNumberGenerators.hpp"

// QUESTION 4 CODE //

double Question4Base::S_T_Euler(double S0, double r, double T, double sigma, double gamma, int m){

    // Assumes/uses an evenly divided time grid of m partitions, S0 > 0 and 0.5 <= gamma <= 1.
    double h = T/m;
    double S_tj = S0;
    for(int i = 0; i < m ; i++){ S_tj = EulerIncrement(S_tj, h, r, sigma, gamma);}

    return S_tj;
}

double Question4Base::S_T_Milstein(double S0,  double r, double T, double sigma, double gamma, int m){

    // Assumes/uses an evenly divided time grid of m partitions, S0 > 0 and 0.5 <= gamma < 1.
    double h = T/m;
    double S_tj = S0;
    for(int i = 0; i < m ; i++){ S_tj = MilsteinIncrement(S_tj, h, r, sigma, gamma);}

    return S_tj;
}

double Question4Base::EulerIncrement(double S_tj, double h, double r, double sigma, double gamma){
    double phi = gamma - 1;
    double drift = ( r - 0.5*sigma*sigma*pow(S_tj,2*phi) )*h;
    double Z = normal_rv(0.0,1.0);
    double diffusion = sigma*pow(S_tj,phi)*sqrt(h)*Z;

    return S_tj*exp(drift + diffusion);
}

double Question4Base::MilsteinIncrement(double S_tj, double h, double r, double sigma, double gamma){
    double phi = gamma - 1;
    double omega = 2*gamma - 3;
    double drift = ( r - 0.5*sigma*sigma*pow(S_tj,2*phi) )*h;
    double Z = normal_rv(0.0,1.0);
    double diffusion = sigma*pow(S_tj,phi)*sqrt(h)*Z;
    double MilsteinTerm = -0.5*sigma*sigma*phi*pow(S_tj,omega)*h*(Z*Z - 1);

    return S_tj*exp(drift + diffusion + MilsteinTerm);
}

CEVMonteCarlo::CEVMonteCarlo(double C, double H, double S0, double r, double T, double sigma, int n, double gamma, int m){
    // Initialising private member variables.
        myC     = C;
        myH     = H;
        myS0    = S0;
        myr     = r;
        myT     = T;
        mysigma = sigma;
        myn     = n;
        mygamma = gamma;
        mym     = m;

}

double CEVMonteCarlo::getEstimate(void){

    double temp, ST_Euler, ST_Milstein, sumEuler = 0, sumEulerSquared = 0, sumMilstein = 0, sumMilsteinSquared = 0;
    for(int i = 0; i < myn; i++){
        // Performing myn = n Monte-Carlo iterations, to produce a numerical estimate of the
        // digital option price.

        // Generating our CEV stock path using the Euler approximation .
        ST_Euler = S_T_Euler(myS0,myr,myT,mysigma,mygamma,mym);

        // Generating our CEV stock path using the Milstein approximation .
        ST_Milstein = S_T_Milstein(myS0,myr,myT,mysigma,mygamma,mym);

        // Our Indicator Functions.
        if( ST_Euler > myH ){
            temp = (myC*exp(-1*myr*myT));
            sumEuler += temp;
            sumEulerSquared += (temp*temp);
        }
        if( ST_Milstein > myH ){
            temp = (myC*exp(-1*myr*myT));
            sumMilstein += temp;
            sumMilsteinSquared += (temp*temp);
        }
    }

    myEulerEstimate = sumEuler/myn;        // Our Euler CEV Monte-Carlo Estimate.
    myMilsteinEstimate = sumMilstein/myn; // Our Euler CEV Monte-Carlo Estimate.


    // Calculating Option asymptotic Variance.
    myEulerVariance = ( sumEulerSquared - (myn*myEulerEstimate*myEulerEstimate) )/(myn*(myn-1));
    myMilsteinVariance = ( sumMilsteinSquared - (myn*myMilsteinEstimate*myMilsteinEstimate) )/(myn*(myn-1)); // assigning value to respective private member variable.

    return myMilsteinEstimate;
}

double CEVMonteCarlo::ConfidenceInterval(double epsilon){
    // Calculating Confidence Interval.
    // Function returns the Milstein interval width.

    double alpha = InverseCumulativeNormal(1.0 - (epsilon/2.0));

    // Calculating Euler Interval
    EulerLeftCI = myEulerEstimate - alpha*sqrt(myEulerVariance);  // assigning values to respective
    EulerRightCI = myEulerEstimate + alpha*sqrt(myEulerVariance); // private member variables.

    // Calculating Milstein Interval
    MilsteinLeftCI = myMilsteinEstimate - alpha*sqrt(myMilsteinVariance);  // assigning values to respective
    MilsteinRightCI = myMilsteinEstimate + alpha*sqrt(myMilsteinVariance); // private member variables.
    EulerCIWidth = EulerRightCI - EulerLeftCI;
    MilsteinCIWidth = MilsteinRightCI - MilsteinLeftCI;


    return MilsteinCIWidth; // The Milstein scheme CI interval width.
}

CEVControlVariates::CEVControlVariates(double C, double H, double S0, double r, double T, double sigma, int n,double gamma, int m){
    // Initialising private member variables.
        myC     = C;
        myH     = H;
        myS0    = S0;
        myr     = r;
        myT     = T;
        mysigma = sigma;
        myn     = n;
        mygamma = gamma;
        mym     = m;

}

double CEVControlVariates::bstarhat(void){

    double numeratorSum = 0, denominatorSum = 0;
    double Xi,Yi;

    // Calculating EX and EY =(approx)= Monte-Carlo(Y).
    double EX = exp(myr*myT)*myS0;
    CEVMonteCarlo CEV(myC,myH,myS0,myr,myT,mysigma,myn,mygamma,mym);
    double EY = CEV.getEstimate();

    for(int i = 0; i < myn; i++){
        Xi = S_T_Milstein(myS0,myr,myT,mysigma,mygamma,mym);
        Yi = 0;
        if( Xi > myH ){
            Yi = myC*exp(-1*myr*myT);
        }
        numeratorSum += ((Xi - EX)*(Yi - EY));
        denominatorSum += ((Xi - EX)*(Xi - EX));
    }
    return numeratorSum/denominatorSum; // Our B*hat estimate.

}

double CEVControlVariates::getEstimate(double b){

    double yib, ST, CVEstimate, sum = 0, squaresum = 0;
    double yi,Xi, EX = exp(myr*myT)*myS0; // E[X] = E[ST].
    for(int i = 0; i < myn; i++){

        // Generating Our stock path.
        ST = S_T_Milstein(myS0,myr,myT,mysigma,mygamma,mym);

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

