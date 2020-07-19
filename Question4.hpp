//
//   Question4.hpp
//  MA323 Project
//
//  Created by Tyler Mitchell on 03/06/2020.
//  Copyright © 2020 Tyler Mitchell. All rights reserved.
//

#ifndef _Question4_hpp
#define _Question4_hpp

#include<iostream>
#include<cmath>
#include<cstdlib>
#include<ctime>

// QUESTION 4 CODE //

// CLASSES //

class Question4Base{

   public:

       virtual double S_T_Euler(double S0, double r, double T, double sigma, double gamma, int m);
       // Computes/Simulates a stock path using the Constant Elasticity of Variance (CEV) model
       // for the price of a risky asset and a first order Euler approximation scheme.

       virtual double S_T_Milstein(double S0, double r, double T, double sigma, double gamma, int m);
       // Computes/Simulates a stock path using the Constant Elasticity of Variance (CEV) model
       // for the price of a risky asset and the second order Milstein approximation scheme.

   private:

        virtual double EulerIncrement(double S_tj, double h, double r, double sigma, double gamma);
        // Calculates the S_t(j+1), i.e. the next step/increment of S_tj, using the Euler Scheme.

        virtual double MilsteinIncrement(double S_tj, double h, double r, double sigma, double gamma);
        // Calculates the S_t(j+1), i.e. the next step/increment of S_tj, using the Milstein Scheme.

};

class CEVMonteCarlo:public Question4Base{

    public:

        CEVMonteCarlo(double C, double H, double S0, double r, double T, double sigma, int n, double gamma, int m);
        // Our Constructor

        virtual double getEstimate(void);
        // Q4.3: computes both the Monte-Carlo estimates of the time-0 price of the cash-or-nothing digital option,
        // using:
        // 1) the CEV model for the risky asset price and the first order Euler approximation.
        // 2) the CEV model for the risky asset price and the second order Milstein approximation.
        // This function also calculate the asymptotic variance of both monte-carlo estimators.
        // Returns the Milstein scheme result.

        double ConfidenceInterval(double epsilon);
        // Q4.3 Calculates our (1-epsilon) asymptotic confidence interval for our both Monte-Carlo estimates
        //-> one using the Euler CEV approximation and the other the Milstein CEV approximation.
        // Note, this returns the width of the (1-epsilon) interval for the Euler scheme. To access the values
        // myEulerVariance, myMilsteinVariance, EulerLeftCI, EulerRightCI, MilsteinLeftCI, MilsteinRightCI and
        // MilsteinCIWidth, which are the digital option's analytical variance and (1-epsilon) confidence interval:
        // [LeftCI,RightCI] and the CI width for the Milstein scheme, run this function and then call these
        // public member variables via:  MC.yourDesiredValue.

        double myEulerVariance, myMilsteinVariance, EulerLeftCI, EulerRightCI;
        double MilsteinLeftCI, MilsteinRightCI, myEulerEstimate, myMilsteinEstimate, EulerCIWidth, MilsteinCIWidth;
        // The public member variables described above.


    private:

        double myC, myH, myS0, myr, myT, mysigma, mygamma;
        int myn, mym;
        // Our private member variables used in calculations.

};

class CEVControlVariates:public Question4Base{

    // All functions utilise the Milstein CEV scheme

    public:

        CEVControlVariates(double C, double H, double S0, double r, double T, double sigma, int n, double gamma, int m);
        // Our Constructor

        double bstarhat(void);
        // Our function for approximating b*

        virtual double getEstimate(double b);
        // Computes the CEV Control Variates estimate of the time-0 price of the cash-or-nothing digital option

        double myvariance;
        // Public member variable which stores the value of the Control Variates estimate
        // of the sample variance of the CV estimator.

        private:

            double myC, myH, myS0, myr, myT, mysigma, mygamma;
            int myn, mym;
            // Our private member variables used in calculations.
};

#include <stdio.h>

#endif /* _Question4_hpp */
