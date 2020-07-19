//
//  Question3.hpp
//  MA323 Project
//
//  Created by Tyler Mitchell on 03/06/2020.
//  Copyright © 2020 Tyler Mitchell. All rights reserved.
//

#ifndef Question3_hpp
#define Question3_hpp


#include<iostream>
#include<cmath>
#include<cstdlib>
#include<ctime>

// QUESTION 3 CODE //

// FUNCTIONS //

double analytical(double C, double H, double S0, double r, double T, double sigma);
// Q3.3a: Computes the t=0 price of the cash-or-nothing digital option using the analytical formula.


// CLASSES //

class Question3Base{

   public:

       virtual double getEstimate(void) = 0;
       // A virtual function common to both derived classes.

       virtual double S_T(double S0, double r, double T, double sigma);
       // Computes/Simulates a stock path modelled upon a Geometric Brownian Motion.

};

class MonteCarlo:public Question3Base{

    public:

        MonteCarlo(double C, double H, double S0, double r, double T, double sigma, int n);
        // Our Constructor

        virtual double getEstimate(void);
        // Q3.3b: computes the Monte-Carlo estimate of the time-0 price of the cash-or-nothing digital option

        double ConfidenceInterval(double epsilon);
        // Q3.3 Calculates our (1-epsilon) asymptotic confidence interval for our Monte-Carlo estimate.
        // Note, this returns the width of the (1-epsilon) interval. To access the values
        // myanalyticalVariance, myleftCI and myrightCI, which are the digital option's analytical
        // variance and (1-epsilon) confidence interval: [myleftCI,myrightCI], run this function
        // and then call these public member variables via:  MC.yourDesiredValue.

        double  myanalyticalVariance, myleftCI, myrightCI;
        // The public member variables described above.


    private:

        double myC, myH, myS0, myr, myT, mysigma, myEstimate;
        int myn;
        // Our private member variables used in calculations.

};

class ControlVariates:public Question3Base{

    public:

        ControlVariates(double C, double H, double S0, double r, double T, double sigma, int n);
        // Our Constructor

        double bstarhat(void);
        // Our function for approximating b*

        virtual double getEstimate(void);
        // Q3.5: computes the Control Variates estimate of the time-0 price of the cash-or-nothing digital option

        double myvariance;
        // Public member variable which stores the value of the Control Variates estimate
        // of the sample variance of the CV estimator.

        private:

            double myC, myH, myS0, myr, myT, mysigma, myEstimate;
            int myn;
            // Our private member variables used in calculations.
};

#include <stdio.h>

#endif /* Question3_hpp */
